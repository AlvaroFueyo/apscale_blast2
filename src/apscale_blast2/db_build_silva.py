"""Database builder: SILVA.

Builds a local BLAST database from SILVA FASTA releases. Taxonomy may be parsed
from headers and/or from an optional external taxonomy table.
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


def load_silva_rank_map(path: str) -> dict[str, str]:
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    mapping = {}
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        for line in handle:
            columns = line.rstrip("\n").split("\t")
            if len(columns) >= 3:
                mapping[columns[0].strip()] = columns[2].strip().lower()
    if not mapping or not any(r == "domain" for r in mapping.values()):
        raise ValueError("Not an official SILVA taxonomy rank map")
    return mapping


def silva_taxonomy_from_headers(fasta_path: str, rank_map_path: str | None = None) -> pd.DataFrame:
    """Map lineage paths to explicit ranks using the matching SILVA release."""
    if not rank_map_path:
        raise ValueError("SILVA header taxonomy requires the matching tax_slv rank map; positional rank guessing is not supported")
    rank_map = load_silva_rank_map(rank_map_path)
    rows = []
    for hdr in _iter_fasta_headers(fasta_path):
        token = hdr.split()[0]
        rest = hdr[len(token) :].strip()
        parts = [p.strip() for p in rest.split(";")] if rest else []
        ranks = dict.fromkeys(["superkingdom", "phylum", "class", "order", "family", "genus", "species"], "")
        for index, value in enumerate(parts):
            lineage = ";".join(parts[:index + 1]) + ";"
            rank = rank_map.get(lineage)
            rank = "superkingdom" if rank == "domain" else rank
            if rank in ranks:
                ranks[rank] = value
        if not ranks["superkingdom"]:
            raise ValueError(f"SILVA header not covered by rank map: {token}")
        # The terminal organism label may supply a species, never an inferred genus.
        from .taxonomy_clean import clean_species
        if parts and clean_species(parts[-1], genus=ranks["genus"]):
            ranks["species"] = parts[-1]
        elif parts and ranks["genus"] and parts[-1].startswith(ranks["genus"] + " "):
            ranks["species"] = parts[-1]
        rows.append([token] + list(ranks.values()))

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
    rank_map_path: str | None = None,
    makeblastdb_exe: str = "makeblastdb",
    keep_source: bool = True,
) -> str:
    """Build a SILVA BLAST database.

    Provide a normalized accession/rank table or an official tax_slv rank map.
    Rank inference by the position of a label is deliberately unsupported.
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
        if taxonomy_path and rank_map_path:
            raise ValueError("Provide taxonomy_path or rank_map_path, not both")
        if rank_map_path:
            tax_df = silva_taxonomy_from_headers(resolved, rank_map_path)
        elif taxonomy_path:
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
            if rank_map_path:
                shutil.copy2(rank_map_path, os.path.join(tmp_out, "source_rank_map" + Path(rank_map_path).suffix))

        # BLAST indices
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
