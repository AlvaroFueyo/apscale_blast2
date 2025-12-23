"""Core execution pipeline.

This module:
- Splits input FASTA files into smaller subsets
- Runs local NCBI BLAST (blastn/megablast) against a selected database
- Applies prefilters and APSCALE-like selection/flag logic
- Writes Excel outputs (raw BLAST hits and taxonomy assignments)

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

BLAST_OUTFMT = "6 qseqid sseqid sacc saccver pident evalue qcovs qcovhsp"

# In-memory cache to avoid reloading taxonomy when multiple FASTA files
# use the same database within a single execution.
# Key: absolute BLAST DB prefix (no extension), normalised.
_TAX_CACHE: Dict[str, Dict[str, List[str]]] = {}

def _tax_cache_key(db_prefix: str) -> str:
    # Normaliza y absolutiza para que rutas equivalentes se cacheen igual.
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
    max_target_seqs: int = 20
    masking: bool = True
    min_qcov: float = 50.0
    prefer_qcov: float = 90.0  # filtro blando por query coverage (se aplica si deja hits)
    max_evalue: float = 1e-3
    min_pident: float = 50.0
    blastn_exe: str = "blastn"
    keep_tsv: bool = False
    thresholds: str = "97,95,90,87,85"
    log_level: str = "INFO"
    inline_perc_identity: bool = True

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
    cols = ["qseqid","sseqid","sacc","saccver","pident","evalue","qcovs","qcovhsp"]
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
        rows.append({
            "unique ID": r["qseqid"],
            "Sequence ID": key or r["sseqid"],
            "Kingdom": (ranks[0] if ranks else ""),
            "Phylum":  (ranks[1] if ranks else ""),
            "Class":   (ranks[2] if ranks else ""),
            "Order":   (ranks[3] if ranks else ""),
            "Family":  (ranks[4] if ranks else ""),
            "Genus":   genus,
            "Species": species,
            "Similarity": float(r["pident"]) if r["pident"]==r["pident"] else 0.0,
            "evalue": float(r["evalue"]) if r["evalue"]==r["evalue"] else 1.0,
            "query_coverage": float(r["qcov"]) if "qcov" in r and r["qcov"]==r["qcov"] else 0.0
        })
    hits = pd.DataFrame(rows, columns=["unique ID","Sequence ID","Kingdom","Phylum","Class","Order","Family","Genus","Species","Similarity","evalue","query_coverage"])

    fasta_dir = os.path.dirname(os.path.abspath(query_fasta)); parent = os.path.dirname(fasta_dir)
    raw_dir = os.path.join(parent, "raw_blast"); tax_dir = os.path.join(parent, "taxonomy")
    os.makedirs(raw_dir, exist_ok=True); os.makedirs(tax_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(query_fasta))[0]
    raw_xlsx = os.path.join(raw_dir, f"{base}_raw_blast.xlsx")
    print("Writing raw output...", flush=True)
    with pd.ExcelWriter(raw_xlsx) as xw:
        hits.to_excel(xw, index=False, sheet_name="raw")

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

        # Filtro "blando" por query coverage: si hay hits con coverage >= prefer_qcov,
        # nos quedamos SOLO con esos. Si no, mantenemos todos los hits del query.
        if float(opts.prefer_qcov) > 0:
            sub_hi = sub[sub["query_coverage"] >= float(opts.prefer_qcov)].copy()
            if not sub_hi.empty:
                sub = sub_hi

        # Hit selection (compatibility with classic APSCALE modes)
        # - modo 1: Similarity -> evalue
        # - modo 2: evalue -> Similarity
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
            row = rr[TAX_COLS + ["Similarity","evalue","query_coverage"]].to_dict()
            trim_by_similarity(row, float(max_sim_ref), thr)
            rows2_full.append({"unique ID": qid, **row})

        rows_species = [r for r in rows2_full if r["Species"]]
        if len(rows_species) == 1:
            r = rows_species[0].copy(); r["Flag"] = ""; r["Ambiguous taxa"] = ""; r["unique ID"] = qid; out_rows.append(r); continue
        if len(rows_species) > 1:
            from collections import Counter
            tax_keys = [(r["Kingdom"], r["Phylum"], r["Class"], r["Order"], r["Family"], r["Genus"], r["Species"]) for r in rows_species]
            counts = Counter(tax_keys)
            maxdup = max(counts.values()) if counts else 0
            if maxdup > 1 and list(counts.values()).count(maxdup) == 1:
                dom_tax = max(counts, key=counts.get)
                dom_row = next(r for r in rows_species
                               if (r["Kingdom"], r["Phylum"], r["Class"], r["Order"], r["Family"], r["Genus"], r["Species"]) == dom_tax)
                dom_row = dom_row.copy()
                dom_row["Flag"] = "F1 (Dominant species)"
                dom_row["Ambiguous taxa"] = ", ".join(sorted(set(r["Species"] for r in rows_species)))
                dom_row["unique ID"] = qid
                out_rows.append(dom_row)
                continue

        seen=set(); rows3=[]
        for r in rows2_full:
            k=(r["Kingdom"],r["Phylum"],r["Class"],r["Order"],r["Family"],r["Genus"],r["Species"])
            if k in seen: continue
            seen.add(k); rows3.append(r)

        amb_species = sorted(set([r["Species"] for r in rows2_full if r["Species"]]))
        chosen = choose_flag_rest(rows3, amb_species)
        chosen["unique ID"] = qid
        out_rows.append(chosen)

    final_df = pd.DataFrame(out_rows, columns=["unique ID"]+TAX_COLS+["Similarity","query_coverage","evalue","Flag","Ambiguous taxa"])
    cat = pd.Categorical(final_df["unique ID"], categories=fasta_order, ordered=True)
    final_df = final_df.assign(_ord=cat).sort_values("_ord", kind="stable").drop(columns=["_ord"])

    out_xlsx = os.path.join(tax_dir, f"{base}_taxonomy.xlsx")
    print("Writing taxonomy output...", flush=True)
    with pd.ExcelWriter(out_xlsx) as xw:
        final_df.to_excel(xw, index=False, sheet_name="Taxonomy table")

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

    print(f"BLAST finalizado en {(time.time()-t0)/60:.1f} min.", flush=True)
    return {"raw_xlsx": raw_xlsx, "filtered_xlsx": out_xlsx}
