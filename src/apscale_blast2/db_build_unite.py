"""Database builder: UNITE.

Builds a local BLAST database from a UNITE FASTA release and extracts taxonomy
from sequence headers.
"""

from __future__ import annotations

import gzip
import os
import shutil
import tempfile
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


def unite_taxonomy_table(fasta_path: str) -> pd.DataFrame:
    """Parse UNITE headers into a taxonomy table.

    Replicates Till's logic:
      taxonomy = record.id.split('|')[4] with k__/p__/...; keys k,p,c,o,f,g,s.
    """
    rows = []
    for hdr in _iter_fasta_headers(fasta_path):
        token = hdr.split()[0]
        parts = token.split("|")
        accession = "|".join(parts[:4])
        tax = parts[4] if len(parts) > 4 else ""
        if not tax or "__" not in tax:
            raise ValueError("Unsupported UNITE header; use the general FASTA release")
        tax_dict = {}
        for item in tax.split(";"):
            if "__" in item:
                k, v = item.split("__", 1)
                tax_dict[k.strip()] = v.replace("_", " ").strip()
        taxonomy = [tax_dict.get(t, "") for t in ["k", "p", "c", "o", "f", "g", "s"]]
        rows.append([accession] + taxonomy)

    return pd.DataFrame(
        rows,
        columns=[
            "Accession",
            "superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
    )


def build_unite_db(
    *,
    input_path: str,
    db_home: str,
    name: str | None = None,
    makeblastdb_exe: str = "makeblastdb",
    keep_source: bool = True,
) -> str:
    """Build an UNITE BLAST database folder and return its path."""
    if not name:
        name = _strip_known_suffixes(Path(input_path).name)
    out_dir = db_folder_for_name(db_home, name)
    if os.path.isdir(out_dir) and os.listdir(out_dir):
        raise FileExistsError(f"The database '{name}' already exists at: {out_dir}")

    tmp_root = tempfile.mkdtemp(prefix="apscale_blast2_build_")
    try:
        tmp_out = os.path.join(tmp_root, "db_build")
        os.makedirs(tmp_out, exist_ok=True)
        tmp_db_dir = os.path.join(tmp_out, "db")
        os.makedirs(tmp_db_dir, exist_ok=True)

        resolved = _resolve_input_file(input_path, tmp_root)
        tax_df = unite_taxonomy_table(resolved)

        if keep_source:
            src_dst = os.path.join(tmp_out, f"source{''.join(Path(resolved).suffixes)}")
            shutil.copy2(resolved, src_dst)

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
