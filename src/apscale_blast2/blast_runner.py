"""Core execution pipeline.

This module:
- Splits input FASTA files into smaller subsets
- Runs local NCBI BLAST (blastn/megablast) against a selected database
- Applies prefilters and APSCALE-like selection/flag logic
- Writes outputs (Excel or Parquet) for raw BLAST hits and taxonomy assignments

Note: BLAST databases are accessed by the external BLAST executables; we only
cache the taxonomy mapping table in-memory within a single Python run.
"""

from __future__ import annotations
import os, subprocess, time, re, sys, logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import List, Dict
import pandas as pd
from tqdm import tqdm

from .io_utils import split_fasta, read_fasta_order
from .dbs import DatabaseSpec, ensure_db_prefix
from .taxmap import load_taxmap_as_dict
from .filtering import thresholds_to_dict, trim_by_similarity, choose_flag_rest
from .taxonomy_clean import clean_species, clean_genus
from . import __version__

# Raw BLAST columns (tabular outfmt 6). Keep this stable unless you bump a major version.
# - sseqid / sacc for easy trace-back to local DB / GenBank
# - mismatch / gapopen for extra QC without bloating exports too much
BLAST_OUTFMT = (
    "6 qseqid sseqid sacc saccver pident evalue qcovs qcovhsp mismatch gapopen"
)

# In-memory cache to avoid reloading taxonomy when multiple FASTA files
# use the same database within a single execution.
# Key: absolute BLAST DB prefix (no extension), normalised.
_TAX_CACHE: Dict[str, Dict[str, List[str]]] = {}

def _tax_cache_key(db_prefix: str) -> str:
    # Normalise and absolutise so that equivalent routes are cached in the same way.
    try:
        p = os.path.abspath(db_prefix)
    except Exception:
        p = db_prefix
    return _norm(p)

def get_tax_dict_cached(db_prefix: str, logger: logging.Logger) -> Dict[str, List[str]]:
    """Load the taxonomy mapping associated with a BLAST DB and cache it in memory."""
    key = _tax_cache_key(db_prefix)
    if key in _TAX_CACHE:
        logger.debug("Taxonomy cache hit for %s", key)
        return _TAX_CACHE[key]
    logger.info("Loading taxonomy...")
    tax_dict = load_taxmap_as_dict(db_prefix)
    _TAX_CACHE[key] = tax_dict
    return tax_dict

@dataclass
class RunOptions:
    threads: int = 0
    # Number of concurrent BLAST workers (subset-level). Default is 1 to avoid oversubscription.
    workers: int = 1
    subset_size: int = 100
    # Default: megablast (fast; suitable for relatively similar amplicons/barcodes).
    task: str = "megablast"
    max_target_seqs: int = 30
    masking: bool = True
    min_qcov: float = 50.0
    prefer_qcov: float = 90.0  # soft query coverage filter (applied only if it leaves at least one hit)
    max_evalue: float = 1e-3
    min_pident: float = 50.0
    blastn_exe: str = "blastn"
    keep_tsv: bool = False
    thresholds: str = "97,95,90,87,85"
    log_level: str = "INFO"
    inline_perc_identity: bool = True

    # Output
    # - excel: write .xlsx using openpyxl (default)
    # - parquet: write .parquet.snappy (requires pyarrow)
    output_format: str = "excel"

    # Flagging / assignment scheme
    # - apscale2: new MRCA-based trimming flags (default)
    # - apscale: legacy APSCALE-BLAST behaviour (dominant-taxon flags, no qcov filters, blastn task)
    flag_scheme: str = "apscale2"

    # Hit selection order *before* applying thresholds/flags.
    # 1 = Similarity -> evalue  ("mode 1" (classic): max similarity first; tie-breaker -> min E-value)
    # 2 = E-value -> Similarity  (min E-value first; tie-breaker -> max similarity)
    filter_mode: int = 1

def _norm(p: str) -> str:
    return p.replace('\\','/')

def _run_single(
    blastn: str,
    db_prefix: str,
    query_fa: str,
    threads: int,
    task: str,
    max_target: int,
    masking: bool,
    log_debug: bool,
    perc_identity: float | None,
    max_evalue: float,
    min_qcov_hsp: float,
) -> str:
    outp = query_fa + ".tsv"
    cmd = [blastn, "-db", _norm(db_prefix), "-query", _norm(query_fa), "-outfmt", BLAST_OUTFMT, "-out", _norm(outp),
           "-task", task, "-max_hsps", "1", "-max_target_seqs", str(max_target), "-num_threads", str(threads)]
    # Early filtering (performance): avoids generating hits that will be discarded later.
    if max_evalue and max_evalue > 0:
        cmd += ["-evalue", str(max_evalue)]
    if min_qcov_hsp and min_qcov_hsp > 0:
        cmd += ["-qcov_hsp_perc", str(min_qcov_hsp)]
    if not masking: cmd += ["-dust","no","-soft_masking","false"]
    if perc_identity is not None and perc_identity > 0:
        cmd += ["-perc_identity", str(perc_identity)]
    if log_debug:
        logging.getLogger("apscale_blast2").debug("CMD: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, capture_output=(not log_debug), text=True)
    except Exception as e:
        logging.getLogger("apscale_blast2").error("Failed to run BLAST: %s", e)
        print("BLAST CMD:", " ".join(cmd), flush=True)
        raise
    return outp

def _before_hash(x: str) -> str:
    return x.split("###",1)[0]

def _before_semicolon(x: str) -> str:
    return x.split(";",1)[0]

def _first_token(x: str) -> str:
    return x.split()[0] if x else x

def _strip_lcl(x: str) -> str:
    return re.sub(r"^lcl\|","", x)

def _strip_version(x: str) -> str:
    return re.sub(r"(\.\d+)$","", x)

def _pipe_core(x: str) -> str:
    parts = x.split("|")
    if len(parts)>=2 and parts[0] in {"gb","ref","emb","sp","tr"}:
        return parts[1]
    return parts[-1] if parts else x

def _resolve_sequence_id(row, tax_dict):
    keys = []
    for k in ["sacc","saccver","sseqid"]:
        v = str(row.get(k,"") or "")
        if v: keys.append(v)
    sseqid = str(row.get("sseqid","") or "")
    sacc   = str(row.get("sacc","") or "")
    saccv  = str(row.get("saccver","") or "")
    keys += [ sseqid, _before_semicolon(sseqid), _before_semicolon(_first_token(sseqid)), _before_hash(sseqid), _first_token(sseqid), _pipe_core(sseqid), _strip_lcl(sseqid),
              _strip_version(sacc), _strip_version(saccv) ]
    for k in keys:
        if k in tax_dict:
            return k, tax_dict[k]
    return None, None


def _extract_accession(subject_id: str) -> str:
    """Best-effort extraction of an accession-like identifier.

    Curated local databases frequently encode taxonomy in the FASTA header
    after a delimiter (commonly "###"), e.g.:

        AP011214.1.70.1027###root_1;Eukaryota_2759;...

    In these cases BLAST may not populate the standard `sacc`/`saccver`
    fields, so we derive a stable identifier from `sseqid`.

    The function:
      1) strips taxonomy suffixes after "###" (if present),
      2) attempts to reduce coordinate-encoded ids to `ACCESSION.VERSION`.

    If it cannot identify an `ACCESSION.VERSION` pattern, it returns the
    stripped id as-is.
    """

    if subject_id is None:
        return ""

    core = str(subject_id).split("###", 1)[0].strip()
    # Many MIDORI/GB-derived headers can look like `AP011214.1.70.1027`.
    # Keep the leading ACCESSION.VERSION if present.
    m = re.match(r"^([A-Za-z]{1,4}\d+\.\d+)", core)
    if m:
        return m.group(1)
    return core


def _try_numeric(series: pd.Series, *, integer: bool) -> pd.Series:
    """Best-effort dtype coercion.

    We try to convert columns that should be numeric (for Parquet exports and
    general downstream sanity). If conversion produces only NaNs while the
    original column had non-empty values, we keep the original dtype to avoid
    breaking edge-case inputs.
    """

    if series is None:
        return series

    s0 = series
    s = pd.to_numeric(s0, errors="coerce")
    # If conversion failed completely but there are non-empty original values,
    # keep as-is.
    try:
        has_nonempty = s0.astype(str).str.strip().ne("").any()
    except Exception:
        has_nonempty = True
    if has_nonempty and s.notna().sum() == 0:
        return s0

    if integer:
        try:
            non_na = s.dropna()
            if not non_na.empty and (non_na % 1 == 0).all():
                return s.round().astype("Int64")
        except Exception:
            pass
    return s.astype(float)


def _prepare_output_dtypes(df: pd.DataFrame, *, int_cols: list[str], float_cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for c in int_cols:
        if c in out.columns:
            out[c] = _try_numeric(out[c], integer=True)
    for c in float_cols:
        if c in out.columns:
            out[c] = _try_numeric(out[c], integer=False)
    return out


def _write_table(df: pd.DataFrame, *, path: str, fmt: str, sheet_name: str | None = None):
    fmt = (fmt or "excel").lower().strip()
    if fmt == "parquet":
        df.to_parquet(path, index=False, compression="snappy")
        return
    # Default: Excel
    with pd.ExcelWriter(path) as xw:
        df.to_excel(xw, index=False, sheet_name=(sheet_name or "Sheet1"))

def _prefilter_hits(raw: pd.DataFrame, min_qcov: float, max_evalue: float, min_pident: float) -> pd.DataFrame:
    df = raw.copy()
    for c in ["pident","evalue","qcovs","qcovhsp"]:
        if c in df.columns: df[c] = pd.to_numeric(df[c], errors="coerce")
    df["qcov"] = df["qcovs"].where(df["qcovs"].notna(), df["qcovhsp"])
    df["qcov"] = pd.to_numeric(df["qcov"], errors="coerce")
    mask = pd.Series(True, index=df.index)
    if min_pident>0: mask &= df["pident"] >= float(min_pident)
    if min_qcov>0:   mask &= df["qcov"]    >= float(min_qcov)
    if max_evalue>0: mask &= df["evalue"]  <= float(max_evalue)
    return df[mask].copy()

def run(query_fasta: str, out_dir: str, db: DatabaseSpec, opts: RunOptions):
    import shutil
    logger = logging.getLogger("apscale_blast2")
    os.makedirs(out_dir, exist_ok=True)
    if opts.threads <= 0:
        try: cpu = os.cpu_count() or 1
        except Exception: cpu = 1
        opts.threads = max(1, cpu - 2)

    fasta_order = read_fasta_order(query_fasta)
    db_prefix = ensure_db_prefix(db.path)
    tax_dict = get_tax_dict_cached(db_prefix, logger)

    subsets_dir = os.path.join(out_dir, "subsets")
    print("Creating FASTA subsets...", flush=True)
    subset_fastas = split_fasta(query_fasta, subsets_dir, subset_size=opts.subset_size)
    total = len(subset_fastas)
    workers = max(1, min(opts.workers, total))
    threads_per = max(1, opts.threads // workers)
    print(f"Starting BLAST: {total} chunks | workers={workers} | total threads={opts.threads}", flush=True)
    print("Running local BLAST (this may take a while)...", flush=True)

    t0 = time.time()
    tsvs: List[str] = []
    pbar = tqdm(total=total, desc=os.path.basename(query_fasta)+" BLAST",
                unit="subset", dynamic_ncols=True, leave=True)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(
                    _run_single,
                    opts.blastn_exe,
                    db_prefix,
                    fa,
                    threads_per,
                    opts.task,
                    opts.max_target_seqs,
                    opts.masking,
                    logger.level<=logging.DEBUG,
                    (opts.min_pident if opts.inline_perc_identity else None),
                    float(opts.max_evalue),
                    float(opts.min_qcov),
                )
                for fa in subset_fastas]
        for fut in as_completed(futs):
            tsvs.append(fut.result())
            pbar.update(1)
    pbar.close()

    print("Starting post-processing...", flush=True)
    print("Merging TSV files...", flush=True)
    cols = [
        "qseqid",
        "sseqid",
        "sacc",
        "saccver",
        "pident",
        "evalue",
        "qcovs",
        "qcovhsp",
        "mismatch",
        "gapopen",
    ]
    dfs = [pd.read_csv(p, sep="\t", names=cols, dtype=str, na_filter=False) for p in tsvs]
    raw = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame(columns=cols)

    print("Applying prefilter...", flush=True)
    raw = _prefilter_hits(raw, opts.min_qcov, opts.max_evalue, opts.min_pident)

    print("Mapping taxonomy...", flush=True)
    rows = []
    for _, r in raw.iterrows():
        key, ranks = _resolve_sequence_id(r, tax_dict)
        if ranks:
            genus = clean_genus(ranks[5])
            species = clean_species(ranks[6])
        else:
            genus = ""; species = ""
        subject_id = key or r["sseqid"]
        # BLAST may not populate sacc/saccver for custom subject ids; also, some
        # curated DB headers embed taxonomy after a delimiter (e.g. "###...").
        # Always normalise to an accession-like identifier for cleaner exports.
        accession_src = r.get("saccver", "") or r.get("sacc", "") or subject_id
        accession = _extract_accession(accession_src)

        # Coverage: prefer qcovs when present; fallback to qcovhsp.
        qcov = None
        try:
            if pd.notna(r.get("qcovs")):
                qcov = float(r["qcovs"])
            elif pd.notna(r.get("qcovhsp")):
                qcov = float(r["qcovhsp"])
        except Exception:
            qcov = None
        if qcov is None or qcov != qcov:
            qcov = 0.0

        rows.append({
            "unique ID": r["qseqid"],
            "Accession": accession,
            "Kingdom": (ranks[0] if ranks else ""),
            "Phylum":  (ranks[1] if ranks else ""),
            "Class":   (ranks[2] if ranks else ""),
            "Order":   (ranks[3] if ranks else ""),
            "Family":  (ranks[4] if ranks else ""),
            "Genus":   genus,
            "Species": species,
            "Similarity": float(r["pident"]) if r["pident"]==r["pident"] else 0.0,
            "evalue": float(r["evalue"]) if r["evalue"]==r["evalue"] else 1.0,
            "query_coverage": float(qcov),
            "mismatch": int(r["mismatch"]) if r.get("mismatch", "") else 0,
            "gapopen": int(r["gapopen"]) if r.get("gapopen", "") else 0,
        })
    hits = pd.DataFrame(
        rows,
        columns=[
            "unique ID",
            "Accession",
            "Kingdom",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species",
            "Similarity",
            "evalue",
            "query_coverage",
            "mismatch",
            "gapopen",
        ],
    )

    # Ensure sane dtypes (especially important for Parquet exports)
    hits = _prepare_output_dtypes(
        hits,
        int_cols=["mismatch", "gapopen"],
        float_cols=["Similarity", "evalue", "query_coverage"],
    )

    fasta_dir = os.path.dirname(os.path.abspath(query_fasta)); parent = os.path.dirname(fasta_dir)
    raw_dir = os.path.join(parent, "raw_blast"); tax_dir = os.path.join(parent, "taxonomy")
    os.makedirs(raw_dir, exist_ok=True); os.makedirs(tax_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(query_fasta))[0]
    fmt = (getattr(opts, "output_format", "excel") or "excel").lower().strip()
    ext = ".parquet.snappy" if fmt == "parquet" else ".xlsx"

    raw_path = os.path.join(raw_dir, f"{base}_raw_blast{ext}")
    print("Writing raw output...", flush=True)
    _write_table(hits, path=raw_path, fmt=fmt, sheet_name="raw")

    print("Applying flags and similarity trimming...", flush=True)
    thr = thresholds_to_dict(opts.thresholds)
    TAX_COLS = ["Kingdom","Phylum","Class","Order","Family","Genus","Species"]

    out_rows: List[Dict[str,object]] = []
    sub_by_q = {qid: sub for qid, sub in hits.groupby("unique ID", sort=False)}
    fasta_ids = read_fasta_order(query_fasta)

    for qid in fasta_ids:
        sub = sub_by_q.get(qid, pd.DataFrame(columns=hits.columns))
        if sub.empty:
            row = {"unique ID": qid}
            for c in TAX_COLS: row[c] = "NoMatch"
            row.update({"Similarity":0.0,"query_coverage":0.0,"evalue":1.0,"Flag":"","Ambiguous taxa":""})
            out_rows.append(row)
            continue

        # Soft query coverage filter: if there are hits with coverage >= prefer_qcov,
        # keep ONLY those; otherwise keep all hits for the query.
        if float(opts.prefer_qcov) > 0:
            sub_hi = sub[sub["query_coverage"] >= float(opts.prefer_qcov)].copy()
            if not sub_hi.empty:
                sub = sub_hi

        # Hit selection (compatibility with classic APSCALE modes)
        # - mode 1: Similarity -> evalue
        # - mode 2: evalue -> Similarity
        if int(opts.filter_mode) == 1:
            max_sim = float(sub["Similarity"].max())
            df1 = sub[(sub["Similarity"] - max_sim).abs() <= 1e-9].copy()
            min_e = float(df1["evalue"].min())
            df1 = df1[df1["evalue"] == min_e].copy()
            max_sim_ref = max_sim
        else:
            min_e = float(sub["evalue"].min())
            df1 = sub[sub["evalue"] == min_e].copy()
            max_sim = float(df1["Similarity"].max()) if not df1.empty else float(sub["Similarity"].max())
            df1 = df1[(df1["Similarity"] - max_sim).abs() <= 1e-9].copy()
            max_sim_ref = max_sim

        rows2_full = []
        for _, rr in df1.iterrows():
            row = rr[TAX_COLS + ["Similarity", "evalue", "query_coverage"]].to_dict()
            trim_by_similarity(row, float(max_sim_ref), thr)
            rows2_full.append({"unique ID": qid, **row})

        # --- Ambiguity handling (F1–F4) ---
        # This follows the APSCALE / apscale_blast logic:
        # 1) Remove duplicate hits after similarity-threshold trimming.
        #    If only one taxon remains, no ambiguity flag is set.
        # 2) If the best similarity is below the species threshold, progressively drop lower
        #    ranks until a single assignment remains (no flag).
        # 3) Otherwise (species-level reference), apply the APSCALE ambiguity flags F1–F4.

        df2 = pd.DataFrame(rows2_full)
        df3 = df2.drop_duplicates().copy()

        # 1) Only one taxon remains -> no ambiguity
        if len(df3) == 1:
            r = df3.iloc[0].to_dict()
            r["Flag"] = ""
            r["Ambiguous taxa"] = ""
            r["unique ID"] = qid
            out_rows.append(r)
            continue

        # 2) Assignment & flags
        if getattr(opts, "flag_scheme", "apscale2") == "apscale":
            # ---- Legacy APSCALE-BLAST flags (dominant-species logic, etc.)

            # Below species threshold -> trim ranks until unique (no flag)
            if float(max_sim_ref) < float(thr["Species"]):
                df_tmp = df3.copy()
                for level in ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom"]:
                    df_tmp[level] = ""
                    df_tmp = df_tmp.drop_duplicates()
                    if len(df_tmp) == 1:
                        break
                r = df_tmp.iloc[0].to_dict()
                r["Flag"] = ""
                r["Ambiguous taxa"] = ""
                r["unique ID"] = qid
                out_rows.append(r)
                continue

            # Species-level reference -> ambiguity flags
            df2c = df2.copy()
            df2c["duplicate_count"] = df2c.groupby(df2c.columns.tolist(), dropna=False).transform("size")
            max_dup = int(df2c["duplicate_count"].max()) if not df2c.empty else 0
            df_dom = df2c[df2c["duplicate_count"] == max_dup].drop_duplicates()

            # F1: dominant species present
            if len(df_dom) == 1:
                dom_row = df_dom.drop(columns=["duplicate_count"]).iloc[0].to_dict()
                dom_row["Flag"] = "F1 (Dominant species)"
                dom_row["Ambiguous taxa"] = ", ".join(sorted({s for s in df2["Species"].drop_duplicates().tolist() if s}))
                dom_row["unique ID"] = qid
                out_rows.append(dom_row)
                continue

            rows3 = df3.to_dict(orient="records")
            amb_species = sorted({s for s in df2["Species"].drop_duplicates().tolist() if s})
            chosen = choose_flag_rest(rows3, amb_species)
            chosen["unique ID"] = qid
            out_rows.append(chosen)
            continue

        # ---- APSCALE-BLAST2 flags (MRCA-based, no dominance)
        # Deduplicate by taxonomy before counting diversity (prevents false flags when equivalent records exist).
        df2_best = df2.copy()
        df2_best["_evalue_num"] = pd.to_numeric(df2_best.get("evalue"), errors="coerce")
        df2_best["_sim_num"] = pd.to_numeric(df2_best.get("Similarity"), errors="coerce")
        # Pick the best hit per unique taxonomy profile
        df2_best = (
            df2_best.sort_values(["_sim_num", "_evalue_num"], ascending=[False, True])
                   .drop_duplicates(subset=TAX_COLS, keep="first")
                   .drop(columns=["_evalue_num", "_sim_num"], errors="ignore")
        )

        # Some reference taxmaps occasionally contain rank placeholders (e.g. "genus")
        # as literal values in taxonomy columns. Treat these as missing.
        _placeholders = {"kingdom", "phylum", "class", "order", "family", "genus", "species"}

        def _clean_tax(v: object) -> str:
            """Normalise taxonomy cell values.

            Drops blanks and rank placeholders such as 'genus' that may appear in
            some reference taxmaps.
            """
            if v is None:
                return ""
            s = str(v).strip()
            if not s:
                return ""
            if s.lower() in _placeholders:
                return ""
            return s

        def _uniq_nonempty(col: str) -> set[str]:
            if col not in df2_best.columns or df2_best.empty:
                return set()
            vals = df2_best[col].fillna("").astype(str).map(_clean_tax)
            return {v for v in vals.tolist() if v}

        ranks_high_to_low = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        ranks_low_to_high = list(reversed(ranks_high_to_low))
        plural = {
            "Species": "species",
            "Genus": "genera",
            "Family": "families",
            "Order": "orders",
            "Class": "classes",
            "Phylum": "phyla",
            "Kingdom": "kingdoms",
        }

        # Find deepest common rank (MRCA) across remaining taxa, ignoring blanks where possible.
        mrca_rank = None
        for rnk in ranks_high_to_low:
            vals = _uniq_nonempty(rnk)
            if len(vals) == 1:
                mrca_rank = rnk
                continue
            # If there are 0 values, taxonomy is incomplete at this rank; stop here.
            if len(vals) == 0:
                break
            # >1 unique values => stop; MRCA is previous (higher) rank.
            break

        # Determine the diversity rank immediately below the MRCA (used for flag numbering/text).
        diversity_rank = None
        if mrca_rank is None:
            diversity_rank = "Kingdom"
        else:
            idx = ranks_high_to_low.index(mrca_rank)
            if idx < len(ranks_high_to_low) - 1:
                diversity_rank = ranks_high_to_low[idx + 1]
            else:
                diversity_rank = None

        # If there is no diversity below MRCA (i.e., only one taxon remains), no flag.
        needs_flag = False
        if diversity_rank is not None:
            vals = _uniq_nonempty(diversity_rank)
            needs_flag = len(vals) > 1

        # Choose a representative/best row for metrics
        best_row = df2_best.copy()
        best_row["_evalue_num"] = pd.to_numeric(best_row.get("evalue"), errors="coerce")
        best_row["_sim_num"] = pd.to_numeric(best_row.get("Similarity"), errors="coerce")
        best_row = best_row.sort_values(["_sim_num", "_evalue_num"], ascending=[False, True]).drop(columns=["_evalue_num", "_sim_num"], errors="ignore")
        out = best_row.iloc[0].to_dict()

        if not needs_flag:
            out["Flag"] = ""
            out["Ambiguous taxa"] = ""
            out["unique ID"] = qid
            out_rows.append(out)
            continue

        # Build flag + ambiguous taxa list.
        #
        # Even if we trim the final assignment up to the MRCA (e.g. to Family),
        # it is useful for manual auditing to keep as much detail as possible
        # about what the surviving hits actually were. Therefore, instead of
        # listing only the values at `diversity_rank` (e.g. just genera), we
        # report each ambiguous taxon at its most specific available rank
        # (Species > Genus > Family > ...).
        amb: list[str] = []
        if diversity_rank:
            for _, rr in df2_best.iterrows():
                genus = _clean_tax(rr.get("Genus", ""))
                species = _clean_tax(rr.get("Species", ""))
                family = _clean_tax(rr.get("Family", ""))
                order = _clean_tax(rr.get("Order", ""))
                klass = _clean_tax(rr.get("Class", ""))
                phylum = _clean_tax(rr.get("Phylum", ""))
                kingdom = _clean_tax(rr.get("Kingdom", ""))

                label = ""
                if species:
                    if genus:
                        # Many taxmaps store Species as a full binomial ("Genus species").
                        if species.lower().startswith((genus + " ").lower()):
                            label = species
                        else:
                            # If Species is just the epithet, prefix the genus; otherwise keep as-is.
                            label = f"{genus} {species}".strip() if " " not in species else species
                    else:
                        label = species
                elif genus:
                    label = genus
                else:
                    label = family or order or klass or phylum or kingdom

                # Guard against placeholders slipping through (best-effort)
                if label and label.lower() not in _placeholders:
                    amb.append(label)
            amb = sorted(set(amb))
        flag_num = ranks_low_to_high.index(diversity_rank) + 1 if diversity_rank else 1
        out["Flag"] = f"Fl{flag_num} Two or more {plural.get(diversity_rank, diversity_rank.lower() + 's')} (trimming to MRCA)"
        out["Ambiguous taxa"] = ", ".join(amb)

        # Trim taxonomy to MRCA (blank all ranks below MRCA). If MRCA unknown, blank everything.
        if mrca_rank is None:
            for c in TAX_COLS:
                out[c] = ""
        else:
            mrca_idx = ranks_high_to_low.index(mrca_rank)
            for c in ranks_high_to_low[mrca_idx + 1:]:
                out[c] = ""

        out["unique ID"] = qid
        out_rows.append(out)

    final_df = pd.DataFrame(out_rows, columns=["unique ID"]+TAX_COLS+["Similarity","query_coverage","evalue","Flag","Ambiguous taxa"])
    cat = pd.Categorical(final_df["unique ID"], categories=fasta_order, ordered=True)
    final_df = final_df.assign(_ord=cat).sort_values("_ord", kind="stable").drop(columns=["_ord"])

    final_df = _prepare_output_dtypes(
        final_df,
        int_cols=[],
        float_cols=["Similarity", "query_coverage", "evalue"],
    )

    out_path = os.path.join(tax_dir, f"{base}_taxonomy{ext}")
    print("Writing taxonomy output...", flush=True)
    _write_table(final_df, path=out_path, fmt=fmt, sheet_name="Taxonomy table")

    # Sidecar run-info for reproducibility (kept small on purpose).
    try:
        from datetime import datetime

        runinfo_path = os.path.join(tax_dir, f"{base}.runinfo.txt")

        # DB timestamp (best-effort)
        db_mtime = "unknown"
        try:
            # Prefer the index file timestamp if we can find it.
            base_dir = os.path.dirname(db_prefix)
            cand = None
            for suf in [".nal", ".nin", ".nhr", ".nsq"]:
                p = db_prefix + suf
                if os.path.exists(p):
                    cand = p
                    break
            cand = cand or (db.path if os.path.exists(db.path) else base_dir)
            db_mtime = datetime.fromtimestamp(os.path.getmtime(cand)).isoformat(timespec="seconds")
        except Exception:
            pass

        lines = [
            f"created_at\t{datetime.now().isoformat(timespec='seconds')}",
            f"tool_version\t{__version__}",
            f"output_format\t{fmt}",
            f"flag_scheme\t{opts.flag_scheme}",
            f"filter_mode\t{getattr(opts, 'filter_mode', 1)}",
            f"task\t{opts.task}",
            f"max_target_seqs\t{opts.max_target_seqs}",
            f"threads\t{opts.threads}",
            f"workers\t{opts.workers}",
            f"subset_size\t{opts.subset_size}",
            f"min_query_coverage\t{opts.min_qcov}",
            f"prefer_query_coverage\t{opts.prefer_qcov}",
            f"max_evalue\t{opts.max_evalue}",
            f"min_pident\t{opts.min_pident}",
            f"thresholds\t{opts.thresholds}",
            f"db_path\t{db.path}",
            f"db_prefix\t{db_prefix}",
            f"db_mtime\t{db_mtime}",
        ]
        with open(runinfo_path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
    except Exception:
        # Never fail a run because of metadata.
        pass

    if not opts.keep_tsv:
        print("Cleaning up temporary files...", flush=True)
        shutil.rmtree(subsets_dir, ignore_errors=True)
        try:
            shutil.rmtree(out_dir, ignore_errors=True)
        except Exception:
            pass
        try:
            root_tmp = os.path.dirname(out_dir)
            if os.path.isdir(root_tmp) and not os.listdir(root_tmp):
                os.rmdir(root_tmp)
        except Exception:
            pass

    print(f"BLAST finished in {(time.time()-t0)/60:.1f} min.", flush=True)
    return {
        "raw_output": raw_path,
        "taxonomy_output": out_path,
        # Backward-compatible keys (paths; may be .xlsx or .parquet.snappy)
        "raw_xlsx": raw_path,
        "filtered_xlsx": out_path,
    }
