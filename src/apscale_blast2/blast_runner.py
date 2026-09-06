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
import os, subprocess, re, logging
from dataclasses import dataclass
import math
from typing import List, Dict
import pandas as pd

from .dbs import DatabaseSpec
from .taxmap import load_taxmap_as_dict
from .filtering import thresholds_to_dict

TAX_COLS = ["Kingdom","Phylum","Class","Order","Family","Genus","Species"]
RANKS_HIGH_TO_LOW = TAX_COLS
RANKS_LOW_TO_HIGH = list(reversed(RANKS_HIGH_TO_LOW))
TAX_PLACEHOLDERS = {"kingdom", "phylum", "class", "order", "family", "genus", "species"}
RANK_PLURALS = {
    "Species": "species",
    "Genus": "genera",
    "Family": "families",
    "Order": "orders",
    "Class": "classes",
    "Phylum": "phyla",
    "Kingdom": "kingdoms",
}

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
_TAX_SIGNATURES = {}

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
    from .taxmap import find_taxmaps_paths
    signature = tuple((p, os.stat(p).st_mtime_ns, os.stat(p).st_size) for p in find_taxmaps_paths(db_prefix))
    if key in _TAX_CACHE and _TAX_SIGNATURES.get(key) == signature:
        logger.debug("Taxonomy cache hit for %s", key)
        return _TAX_CACHE[key]
    logger.info("Loading taxonomy...")
    tax_dict = load_taxmap_as_dict(db_prefix)
    if len(_TAX_CACHE) >= 4:
        _TAX_CACHE.clear()
        _TAX_SIGNATURES.clear()
    _TAX_CACHE[key] = tax_dict
    _TAX_SIGNATURES[key] = signature
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
    prefer_qcov: float = 75.0  # soft query coverage filter (applied only if it leaves at least one hit)
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
    output_dir: str | None = None
    overwrite: bool = False

    # Flagging / assignment scheme
    # - apscale2: new MRCA-based trimming flags (default)
    # - apscale: legacy APSCALE-BLAST behaviour (dominant-taxon flags, no qcov filters, blastn task)
    flag_scheme: str = "apscale2"

    # Hit selection order *before* applying thresholds/flags.
    # 1 = Similarity -> evalue  ("mode 1" (classic): max similarity first; tie-breaker -> min E-value)
    # 2 = E-value -> Similarity  (min E-value first; tie-breaker -> max similarity)
    filter_mode: int = 1

    def __post_init__(self):
        thresholds_to_dict(self.thresholds)
        for name, lo, hi in [("threads", 0, 1024), ("workers", 1, 128), ("subset_size", 1, 1000000), ("max_target_seqs", 1, 1000)]:
            value = getattr(self, name)
            if not isinstance(value, int) or not lo <= value <= hi:
                raise ValueError(f"{name} must be an integer in [{lo}, {hi}]")
        for name, lo, hi in [("min_qcov", 0, 100), ("prefer_qcov", 0, 100), ("min_pident", 0, 100), ("max_evalue", 1e-300, 1)]:
            value = getattr(self, name)
            if not math.isfinite(value) or not lo <= value <= hi:
                raise ValueError(f"Invalid {name}: {value}")
        if self.flag_scheme not in {"apscale", "apscale2"} or self.filter_mode not in {1, 2}:
            raise ValueError("Invalid assignment scheme or filter mode")
        if self.task not in {"blastn", "megablast"} or self.output_format not in {"excel", "parquet"}:
            raise ValueError("Invalid BLAST task or output format")
        if self.flag_scheme == "apscale":
            self.task, self.min_qcov, self.prefer_qcov, self.max_target_seqs = "blastn", 0.0, 0.0, 20

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
    cancel_event=None,
) -> str:
    if cancel_event is not None and cancel_event.is_set():
        raise InterruptedError("BLAST cancelled")
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
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, errors="replace") as process:
            try:
                while True:
                    if cancel_event is not None and cancel_event.is_set():
                        raise InterruptedError("BLAST cancelled")
                    try:
                        stdout, stderr = process.communicate(timeout=0.25)
                        break
                    except subprocess.TimeoutExpired:
                        continue
                if process.returncode:
                    raise subprocess.CalledProcessError(process.returncode, cmd, output=stdout, stderr=stderr)
                if log_debug and stderr:
                    logging.getLogger("apscale_blast2").debug("BLAST: %s", stderr.strip())
            except BaseException:
                process.terminate()
                try:
                    process.communicate(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.communicate()
                raise
    except Exception as e:
        logging.getLogger("apscale_blast2").error("Failed to run BLAST: %s; stderr: %s", e, getattr(e, "stderr", ""))
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
    keys += [ sseqid, sseqid.split("|k__", 1)[0], _before_semicolon(sseqid), _before_semicolon(_first_token(sseqid)), _before_hash(sseqid), _first_token(sseqid), _pipe_core(sseqid), _strip_lcl(sseqid),
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
    """Run one FASTA; out_dir is retained for API compatibility, never deleted.

    Outputs follow opts.output_dir or the historical parent-of-FASTA layout.
    Intermediate files live in a uniquely owned output-side directory.
    """
    from .streaming import run as run_streaming
    return run_streaming(query_fasta, out_dir, db, opts)
