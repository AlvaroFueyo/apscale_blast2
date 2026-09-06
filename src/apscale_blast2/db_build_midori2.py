"""Database builder: MIDORI2.

Builds a local BLAST database from a MIDORI2 FASTA (optionally compressed or
packaged). Produces BLAST indices plus a taxonomy mapping table.
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


def _is_gz(p: str) -> bool:
    return p.lower().endswith(".gz")


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


def _midori_token_to_name(t: str) -> str:
    parts = [p for p in str(t).split("_") if p]
    if not parts:
        return ""
    if len(parts) == 1:
        return parts[0]

    # If the last token is a taxid, drop it.
    if parts[-1].isdigit():
        parts = parts[:-1]

    if not parts:
        return ""

    # Higher ranks often stay in one token; species commonly need the first two.
    if len(parts) == 1:
        return parts[0]
    return " ".join(parts)


def _iter_fasta_headers(path: str):
    opener = gzip.open if _is_gz(path) else open
    mode = "rt" if _is_gz(path) else "r"
    with opener(path, mode, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].strip()


def midori2_taxonomy_table(fasta_path: str) -> pd.DataFrame:
    """Parse MIDORI2 headers into a taxonomy table.

    MIDORI2 BLAST+ FASTA headers include taxonomy separated by semicolons.
    This replicates the parsing logic used in Till Macher's builder scripts.
    """
    rows = []
    for hdr in _iter_fasta_headers(fasta_path):
        # Accession = first token before whitespace. Taxonomy follows after ';'.
        tokens = hdr.split()
        if not tokens:
            continue
        first_token = tokens[0]
        parts = first_token.split(';')
        if not parts:
            continue
        accession = parts[0].lstrip('>')
        tax_parts = parts[1:]
        if "###root_" not in accession or len(tax_parts) != 7:
            raise ValueError("Unsupported MIDORI2 header; use the BLAST-formatted FASTA with root plus seven ranks")

        taxonomy = []
        for t in tax_parts:
            taxonomy.append(_midori_token_to_name(t))

        # Ensure exactly 7 ranks (superkingdom..species)
        taxonomy = (taxonomy + [''] * 7)[:7]
        rows.append([accession] + taxonomy)

    return pd.DataFrame(
        rows,
        columns=[
            'Accession',
            'superkingdom',
            'phylum',
            'class',
            'order',
            'family',
            'genus',
            'species',
        ],
    )


def build_midori2_db(
    *,
    input_path: str,
    db_home: str,
    name: str | None = None,
    makeblastdb_exe: str = "makeblastdb",
    keep_source: bool = True,
) -> str:
    """Build a MIDORI2 BLAST database folder and return its path.

    Output layout:
      <db_home>/db_<name>/db/ (BLAST indices)
      <db_home>/db_<name>/db_taxonomy.parquet.snappy
    """
    if not name:
        name = _strip_known_suffixes(Path(input_path).name)
    out_dir = db_folder_for_name(db_home, name)

    # atomic build into temp and then move
    if os.path.isdir(out_dir) and os.listdir(out_dir):
        raise FileExistsError(f"The database '{name}' already exists at: {out_dir}")

    tmp_root = tempfile.mkdtemp(prefix="apscale_blast2_build_")
    try:
        tmp_out = os.path.join(tmp_root, "db_build")
        os.makedirs(tmp_out, exist_ok=True)
        tmp_db_dir = os.path.join(tmp_out, "db")
        os.makedirs(tmp_db_dir, exist_ok=True)

        # Resolve input (handles ZIP)
        resolved = _resolve_input_file(input_path, tmp_root)

        # Parse taxonomy
        tax_df = midori2_taxonomy_table(resolved)

        # Optionally store source file in output
        if keep_source:
            src_dst = os.path.join(tmp_out, f"source{''.join(Path(resolved).suffixes)}")
            shutil.copy2(resolved, src_dst)

        # Build BLAST indices
        prefix = os.path.join(tmp_db_dir, "db")
        from .db_build_common import run_makeblastdb
        run_makeblastdb(makeblastdb_exe, resolved, prefix, tax_df)

        # Write taxonomy parquet
        tax_path = os.path.join(tmp_out, "db_taxonomy.parquet.snappy")
        tax_df.to_parquet(tax_path)

        # Move atomically
        os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        from .db_build_common import install_built_database
        return install_built_database(tmp_out, out_dir)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
