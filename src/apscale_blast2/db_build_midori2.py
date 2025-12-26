"""Database builder: MIDORI2.

Builds a local BLAST database from a MIDORI2 FASTA (optionally compressed or
packaged). Produces BLAST indices plus a taxonomy mapping table.
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
    """Return a local path to a fasta(.gz) file.

    If input is a ZIP, extracts the first fasta-like file.
    """
    p = os.path.abspath(os.path.expanduser(input_path))
    if not os.path.exists(p):
        raise FileNotFoundError(p)
    if p.lower().endswith(".zip"):
        with zipfile.ZipFile(p, "r") as zf:
            members = [m for m in zf.namelist() if not m.endswith("/")]
            # Extract a FASTA-like member from the archive. If no FASTA is present, this is
            # likely a precompiled BLAST database archive (indices only) rather than a source FASTA.
            def is_fasta_member(m: str) -> bool:
                ml = m.lower()
                return any(ml.endswith(ext) for ext in FA_EXTS) or any(ml.endswith(ext + ".gz") for ext in FA_EXTS)

            fasta_members = [m for m in members if is_fasta_member(m)]
            if not fasta_members:
                raise ValueError(
                    "No FASTA file was found inside the ZIP archive. If you intended to install a precompiled "
                    "BLAST database (e.g., a db_*.zip bundle with .nsq/.nin files), use the precompiled database "
                    "installer instead of the MIDORI2 builder."
                )

            # Prefer compressed FASTA members first, then plain FASTA.
            def score(m: str) -> int:
                ml = m.lower()
                if any(ml.endswith(ext + ".gz") for ext in FA_EXTS):
                    return 0
                return 1

            fasta_members = sorted(fasta_members, key=score)
            chosen = fasta_members[0]
            out = os.path.join(workdir, os.path.basename(chosen))
            zf.extract(chosen, workdir)
            # if nested folders, move to out
            extracted = os.path.join(workdir, chosen)
            if extracted != out:
                os.makedirs(os.path.dirname(out), exist_ok=True)
                shutil.move(extracted, out)
                # cleanup empty dirs
                try:
                    root = os.path.join(workdir, os.path.dirname(chosen))
                    if root and os.path.isdir(root):
                        shutil.rmtree(root, ignore_errors=True)
                except Exception:
                    pass
            return out
    return p


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

        taxonomy = []
        for t in tax_parts:
            record_split = t.split('_')
            if len(record_split) == 2:
                taxonomy.append(record_split[0])
            else:
                taxonomy.append(' '.join(record_split[:2]))

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
        cmd = [makeblastdb_exe, "-in", resolved, "-title", "db", "-dbtype", "nucl", "-out", prefix]
        subprocess.run(cmd, check=True)

        # Write taxonomy parquet
        tax_path = os.path.join(tmp_out, "db_taxonomy.parquet.snappy")
        tax_df.to_parquet(tax_path)

        # Move atomically
        os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        shutil.move(tmp_out, out_dir)
        return out_dir
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
