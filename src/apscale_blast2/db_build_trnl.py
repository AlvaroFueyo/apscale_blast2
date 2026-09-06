"""Database builder: trnL.

Builds a local BLAST database for trnL workflows from a FASTA plus an external
taxonomy table.
"""

from __future__ import annotations

import gzip
import os
import shutil
import tempfile
import csv
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
    from .db_build_common import resolve_input_file
    return resolve_input_file(input_path, workdir)


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
        df = pd.read_excel(p, dtype=str)
    else:
        with open(p, encoding="utf-8-sig", newline="") as handle:
            sample = handle.read(8192)
            handle.seek(0)
            if not sample.strip():
                raise ValueError("Empty trnL taxonomy table")
            delimiter = "\t" if "\t" in sample.splitlines()[0] else csv.Sniffer().sniff(sample, delimiters=",;").delimiter
            rows = list(csv.reader(handle, delimiter=delimiter))
        if rows and len(rows[0]) == 2:
            # Official CRUX: no header, accession TAB seven semicolon-separated ranks.
            rows = [[row[0], *row[1].split(";")] for row in rows if row]
        columns = ["Accession", "superkingdom", "phylum", "class", "order", "family", "genus", "species"]
        if not rows:
            raise ValueError("Empty trnL taxonomy table")
        header = [name.strip().lower() for name in rows[0]]
        if header[0] in {"accession", "sequence id", "sequence_id", "id"}:
            df = pd.DataFrame(rows[1:], columns=header)
            df = df.rename(columns={header[0]: "Accession"})
            if not set(columns).issubset(df.columns):
                raise ValueError("trnL table requires named Accession and seven rank columns")
            df = df[columns]
        else:
            if any(len(row) != 8 for row in rows):
                raise ValueError("Headerless CRUX taxonomy must have exactly eight fields")
            df = pd.DataFrame(rows, columns=columns)
    df = df.fillna("")
    columns = [
        "Accession",
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    renamed = {c: str(c).strip().lower() for c in df.columns}
    renamed = {c: "Accession" if name in {"accession", "sequence id", "sequence_id", "id"} else "superkingdom" if name in {"kingdom", "domain"} else name for c, name in renamed.items()}
    df = df.rename(columns=renamed)
    if not set(columns).issubset(df.columns):
        raise ValueError("trnL taxonomy requires named Accession and seven rank columns")
    return df[columns].copy()


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
        from .db_build_common import run_makeblastdb
        run_makeblastdb(makeblastdb_exe, resolved, prefix, tax_df)

        tax_path = os.path.join(tmp_out, "db_taxonomy.parquet.snappy")
        tax_df.to_parquet(tax_path)

        os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        from .db_build_common import install_built_database
        return install_built_database(tmp_out, out_dir)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
