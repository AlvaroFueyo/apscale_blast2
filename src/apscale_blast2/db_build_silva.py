"""Database builder: SILVA.

Builds a local BLAST database from SILVA FASTA releases. Taxonomy may be parsed
from headers and/or from an optional external taxonomy table.
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
    """Return a local path to a fasta(.gz) file. If input is a ZIP, extracts the first fasta-like file."""
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
    """Load a taxonomy table and normalise to Accession+7 ranks if possible."""
    p = os.path.abspath(os.path.expanduser(taxonomy_path))
    if not os.path.exists(p):
        raise FileNotFoundError(p)

    low = p.lower()
    if low.endswith(".parquet") or ".parquet" in os.path.basename(low):
        df = pd.read_parquet(p)
    elif low.endswith(".xlsx") or low.endswith(".xls"):
        df = pd.read_excel(p)
    elif low.endswith(".tsv"):
        df = pd.read_csv(p, sep="\t", dtype=str)
    else:
        df = pd.read_csv(p, dtype=str)

    df = df.fillna("")
    lowcols = {c.lower(): c for c in df.columns}

    # ID column
    id_col = None
    for cand in ["Accession", "Sequence ID", "sequence_id", "id", "header"]:
        if cand in df.columns:
            id_col = cand
            break
        if cand.lower() in lowcols:
            id_col = lowcols[cand.lower()]
            break
    if id_col is None:
        raise ValueError("No ID column was found in the taxonomy table (Accession/Sequence ID).")

    def pick(*names: str) -> str | None:
        for nm in names:
            if nm in df.columns:
                return nm
            if nm.lower() in lowcols:
                return lowcols[nm.lower()]
        return None

    colmap = {
        "superkingdom": pick("superkingdom", "kingdom", "domain", "kingdom_name", "superkingdom_name"),
        "phylum": pick("phylum", "division", "phylum_name"),
        "class": pick("class", "class_name"),
        "order": pick("order", "order_name"),
        "family": pick("family", "family_name"),
        "genus": pick("genus", "genus_name"),
        "species": pick("species", "species_name"),
    }

    out = pd.DataFrame({"Accession": df[id_col].astype(str)})
    for k, c in colmap.items():
        out[k] = df[c].astype(str) if c else ""
    return out[["Accession", "superkingdom", "phylum", "class", "order", "family", "genus", "species"]]


def silva_taxonomy_from_headers(fasta_path: str) -> pd.DataFrame:
    """Attempt to parse SILVA taxonomy from FASTA header lines.

    Many SILVA exports include a semicolon-separated lineage after the first token.
    We map to the APSCALE ranks using a simple heuristic.
    """
    rows = []
    for hdr in _iter_fasta_headers(fasta_path):
        token = hdr.split()[0]
        rest = hdr[len(token) :].strip()
        parts = [p.strip() for p in rest.split(";") if p.strip()] if rest else []

        if len(parts) >= 2:
            superkingdom = parts[0] if len(parts) > 0 else ""
            phylum = parts[1] if len(parts) > 1 else ""
            clazz = parts[2] if len(parts) > 2 else ""
            order = parts[3] if len(parts) > 3 else ""
            family = parts[4] if len(parts) > 4 else ""
            genus = parts[-2] if len(parts) >= 2 else ""
            species = parts[-1] if len(parts) >= 1 else ""
        else:
            superkingdom = phylum = clazz = order = family = genus = species = ""

        rows.append([token, superkingdom, phylum, clazz, order, family, genus, species])

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


def build_silva_db(
    *,
    input_path: str,
    db_home: str,
    name: str | None = None,
    taxonomy_path: str | None = None,
    makeblastdb_exe: str = "makeblastdb",
    keep_source: bool = True,
) -> str:
    """Build a SILVA BLAST database.

    If taxonomy_path is provided, uses it as mapping table. Otherwise tries to parse
    taxonomy from FASTA headers.
    """
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

        # Build taxonomy table
        if taxonomy_path:
            tax_df = _load_taxonomy_table(taxonomy_path)
            # Reduce to sequences present in FASTA, if possible
            try:
                accessions = {hdr.split()[0] for hdr in _iter_fasta_headers(resolved)}
                if accessions:
                    tax_df = tax_df[tax_df["Accession"].astype(str).isin(accessions)].copy()
            except Exception:
                pass
        else:
            tax_df = silva_taxonomy_from_headers(resolved)

        if keep_source:
            src_dst = os.path.join(tmp_out, f"source{''.join(Path(resolved).suffixes)}")
            shutil.copy2(resolved, src_dst)
            if taxonomy_path and os.path.exists(taxonomy_path):
                shutil.copy2(os.path.abspath(os.path.expanduser(taxonomy_path)), os.path.join(tmp_out, "source_taxonomy" + Path(taxonomy_path).suffix))

        # BLAST indices
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
