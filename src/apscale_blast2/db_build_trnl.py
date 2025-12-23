"""Database builder: trnL.

Builds a local BLAST database for trnL workflows from a FASTA plus an external
taxonomy table.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

import pandas as pd

from .db_home import db_folder_for_name


FA_EXTS = (".fa", ".fasta", ".fna")


def _strip_known_suffixes(name: str) -> str:
    n = name
    for suf in [".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna", ".gz", ".zip"]:
        if n.lower().endswith(suf):
            n = n[: -len(suf)]
            break
    return n


def _resolve_input_file(input_path: str, workdir: str) -> str:
    """Return a local path to a fasta(.gz) file.

    If input is a ZIP, extracts the first fasta-like file found.
    """
    p = os.path.abspath(os.path.expanduser(input_path))
    if not os.path.exists(p):
        raise FileNotFoundError(p)
    if p.lower().endswith(".zip"):
        with zipfile.ZipFile(p, "r") as zf:
            members = [m for m in zf.namelist() if not m.endswith("/")]

            def score(m: str) -> int:
                ml = m.lower()
                if any(ml.endswith(ext + ".gz") for ext in FA_EXTS):
                    return 0
                if any(ml.endswith(ext) for ext in FA_EXTS):
                    return 1
                return 2

            members = sorted(members, key=score)
            if not members:
                raise ValueError("The ZIP archive is empty.")
            chosen = members[0]
            zf.extract(chosen, workdir)
            extracted = os.path.join(workdir, chosen)
            out = os.path.join(workdir, os.path.basename(chosen))
            if extracted != out:
                os.makedirs(os.path.dirname(out), exist_ok=True)
                shutil.move(extracted, out)
                try:
                    root = os.path.join(workdir, os.path.dirname(chosen))
                    if root and os.path.isdir(root):
                        shutil.rmtree(root, ignore_errors=True)
                except Exception:
                    pass
            return out
    return p


def _iter_fasta_headers(path: str):
    is_gz = path.lower().endswith(".gz")
    opener = gzip.open if is_gz else open
    mode = "rt" if is_gz else "r"
    with opener(path, mode, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].strip()


def _load_taxonomy_table(taxonomy_path: str) -> pd.DataFrame:
    p = os.path.abspath(os.path.expanduser(taxonomy_path))
    if not os.path.exists(p):
        raise FileNotFoundError(p)

    low = p.lower()
    if low.endswith(".xlsx") or low.endswith(".xls"):
        df = pd.read_excel(p)
    else:
        # Till: sep='[;\t]' engine='python'
        df = pd.read_csv(p, sep=r"[;\t]", engine="python", dtype=str)
    df = df.fillna("")
    if df.shape[1] < 8:
        raise ValueError("The trnL taxonomy table must have at least 8 columns (Accession + 7 ranks).")
    df = df.iloc[:, :8].copy()
    df.columns = [
        "Accession",
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    return df


def build_trnl_db(
    *,
    input_fasta: str,
    taxonomy_path: str,
    db_home: str,
    name: str | None = None,
    makeblastdb_exe: str = "makeblastdb",
    keep_source: bool = True,
) -> str:
    """Build a trnL BLAST database folder and return its path."""
    if not name:
        name = _strip_known_suffixes(Path(input_fasta).name)
    out_dir = db_folder_for_name(db_home, name)
    if os.path.isdir(out_dir) and os.listdir(out_dir):
        raise FileExistsError(f"The database '{name}' already exists at: {out_dir}")

    tmp_root = tempfile.mkdtemp(prefix="apscale_blast2_build_")
    try:
        tmp_out = os.path.join(tmp_root, "db_build")
        os.makedirs(tmp_out, exist_ok=True)
        tmp_db_dir = os.path.join(tmp_out, "db")
        os.makedirs(tmp_db_dir, exist_ok=True)

        resolved = _resolve_input_file(input_fasta, tmp_root)
        tax_df = _load_taxonomy_table(taxonomy_path)

        # Reduce size: keep only accessions present in the FASTA (when possible)
        try:
            accessions = {hdr.split()[0] for hdr in _iter_fasta_headers(resolved)}
            if accessions:
                tax_df = tax_df[tax_df["Accession"].astype(str).isin(accessions)].copy()
        except Exception:
            pass

        if keep_source:
            src_dst = os.path.join(tmp_out, f"source{''.join(Path(resolved).suffixes)}")
            shutil.copy2(resolved, src_dst)
            shutil.copy2(os.path.abspath(os.path.expanduser(taxonomy_path)), os.path.join(tmp_out, "source_taxonomy" + Path(taxonomy_path).suffix))

        prefix = os.path.join(tmp_db_dir, "db")
        cmd = [makeblastdb_exe, "-in", resolved, "-title", "db", "-dbtype", "nucl", "-out", prefix]
        subprocess.run(cmd, check=True)

        tax_path = os.path.join(tmp_out, "db_taxonomy.parquet.snappy")
        tax_df.to_parquet(tax_path)

        os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        shutil.move(tmp_out, out_dir)
        return out_dir
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
