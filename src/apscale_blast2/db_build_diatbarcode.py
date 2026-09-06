"""Database builder: DiatBarcode.

Builds a local BLAST database from the DiatBarcode Excel release.
"""

from __future__ import annotations

import os
import json
import logging
import shutil
import tempfile
from pathlib import Path

import pandas as pd

from .db_home import db_folder_for_name


RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]


def read_diatbarcode(path: str, classification: str = "RCM") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read legacy flat releases or the sequence/tree layout introduced in v16."""
    if classification not in {"RCM", "Kociolek"}:
        raise ValueError("DiatBarcode classification must be RCM or Kociolek")
    with pd.ExcelFile(path) as book:
        if "sequences_info" in book.sheet_names:
            sequences = pd.read_excel(book, sheet_name="sequences_info", dtype=str).fillna("")
            tree = pd.read_excel(book, sheet_name=f"taxo_{classification}", dtype=str).fillna("")
            nodes = {}
            for name, rank, parent in tree[["taxon name", "rank", "parent taxon name"]].itertuples(index=False, name=None):
                name, rank, parent = name.strip(), rank.strip().lower(), parent.strip()
                if not name:
                    continue
                nodes.setdefault(name, set()).add((rank, parent))

            cache = {}

            def paths(node, visited=frozenset()):
                if not node:
                    return [{**dict.fromkeys(RANKS, ""), "missing_nodes": ""}]
                if node in visited:
                    raise ValueError(f"Cycle in DiatBarcode taxonomy: {node}")
                if node in cache:
                    return cache[node]
                if node not in nodes:
                    return [{**dict.fromkeys(RANKS, ""), "missing_nodes": node}]
                result = []
                for rank, parent in sorted(nodes[node]):
                    for ancestor in paths(parent, visited | {node}):
                        lineage = ancestor.copy()
                        if rank in {"domain", "kingdom", "subkingdom"}:
                            lineage["superkingdom"] = node
                        elif rank in lineage:
                            lineage[rank] = node
                        result.append(lineage)
                cache[node] = result
                return result

            rows = []
            for acc, species in sequences[["Sequence ID", "Species"]].itertuples(index=False, name=None):
                alternatives = paths(species.strip())
                consensus, conflicts = dict.fromkeys(RANKS, ""), []
                for rank in RANKS:
                    values = {lineage[rank] for lineage in alternatives if lineage[rank]}
                    if len(values) > 1:
                        conflicts = [f"{rank}: {' / '.join(sorted(values))}"]
                        break
                    if values:
                        consensus[rank] = next(iter(values))
                rows.append({"Accession": acc.strip(), **consensus, "species_original": species.strip(),
                             "taxonomy_conflict": "; ".join(conflicts),
                             "taxonomy_missing_nodes": "; ".join(sorted({r["missing_nodes"] for r in alternatives} - {""})),
                             "source_lineages": json.dumps(alternatives, sort_keys=True) if conflicts or any(r["missing_nodes"] for r in alternatives) else ""})
        else:
            if classification != "RCM":
                raise ValueError("Legacy DiatBarcode workbooks only provide RCM taxonomy")
            sheet = "diatbarcode v12" if "diatbarcode v12" in book.sheet_names else book.sheet_names[0]
            sequences = pd.read_excel(book, sheet_name=sheet, dtype=str).fillna("")
            columns = ["Sequence ID", "Subkingdom (following Algaebase 2018)",
                       "Phylum (following Algaebase 2018)", "Class (following Round, Crawford & Mann 1990)",
                       "Order (following Round, Crawford & Mann 1990)", "Family (following Round, Crawford & Mann 1990)",
                       "Genus", "Species"]
            rows = [dict(zip(["Accession"] + RANKS, row)) for row in sequences[columns].values.tolist()]
    if not {"Sequence ID", "Sequence"}.issubset(sequences.columns):
        raise ValueError("DiatBarcode requires Sequence ID and Sequence columns")
    sequences = sequences.copy()
    sequences["gaps_removed"] = sequences["Sequence"].str.count("-")
    sequences["Sequence"] = sequences["Sequence"].str.replace(r"\s+|-", "", regex=True)
    taxonomy = pd.DataFrame(rows)
    for column in ["taxonomy_conflict", "taxonomy_missing_nodes", "source_lineages"]:
        if column not in taxonomy:
            taxonomy[column] = ""
    if "species_original" not in taxonomy:
        taxonomy["species_original"] = taxonomy["species"]
    conflicts = taxonomy["taxonomy_conflict"].ne("").sum()
    if conflicts:
        logging.getLogger("apscale_blast2").warning("DiatBarcode %s: %s references trimmed to compatible ranks; see build_audit.csv", classification, conflicts)
    missing = taxonomy["taxonomy_missing_nodes"].ne("").sum()
    if missing:
        logging.getLogger("apscale_blast2").warning("DiatBarcode %s: %s references have missing tree ancestors; known ranks retained, see build_audit.csv", classification, missing)
    return sequences, taxonomy


def build_diatbarcode_db(
    *,
    xlsx_path: str,
    db_home: str,
    name: str | None = None,
    makeblastdb_exe: str = "makeblastdb",
    keep_source: bool = True,
    classification: str = "RCM",
) -> str:
    """Build a DiatBarcode BLAST database from the official XLSX release.

    Supports the legacy flat table and v16 tree sheets. RCM is the default;
    Kociolek is selectable only when the workbook contains that classification.
    """
    p = os.path.abspath(os.path.expanduser(xlsx_path))
    if not os.path.exists(p):
        raise FileNotFoundError(p)

    if not name:
        base = Path(p).stem.replace(" ", "_")
        name = f"diatbarcode_{base}" if not base.lower().startswith("diat") else base

    out_dir = db_folder_for_name(db_home, name)
    if os.path.isdir(out_dir) and os.listdir(out_dir):
        raise FileExistsError(f"The database '{name}' already exists at: {out_dir}")

    tmp_root = tempfile.mkdtemp(prefix="apscale_blast2_build_")
    try:
        tmp_out = os.path.join(tmp_root, "db_build")
        os.makedirs(tmp_out, exist_ok=True)
        tmp_db_dir = os.path.join(tmp_out, "db")
        os.makedirs(tmp_db_dir, exist_ok=True)

        diat_df, tax_df = read_diatbarcode(p, classification)
        audit = tax_df[["Accession", "species_original", "taxonomy_conflict", "taxonomy_missing_nodes", "source_lineages"]].copy()
        audit["classification"] = classification
        audit["gaps_removed"] = diat_df["gaps_removed"].to_numpy()
        audit.to_csv(os.path.join(tmp_out, "build_audit.csv"), index=False)

        # FASTA
        fasta_path = os.path.join(tmp_root, "diatbarcode.fasta")
        with open(fasta_path, "w", encoding="utf-8", errors="replace") as fh:
            if "Sequence ID" not in diat_df.columns or "Sequence" not in diat_df.columns:
                raise ValueError("El XLSX de DiatBarcode debe contener columnas 'Sequence ID' y 'Sequence'.")
            for seq_id, seq in diat_df[["Sequence ID", "Sequence"]].values.tolist():
                seq_id = str(seq_id).strip()
                seq = str(seq).strip()
                if not seq_id:
                    continue
                fh.write(f">{seq_id}\n")
                fh.write(f"{seq}\n")

        if keep_source:
            shutil.copy2(p, os.path.join(tmp_out, "source.xlsx"))

        # BLAST indices
        prefix = os.path.join(tmp_db_dir, "db")
        from .db_build_common import run_makeblastdb
        run_makeblastdb(makeblastdb_exe, fasta_path, prefix, tax_df)

        tax_path = os.path.join(tmp_out, "db_taxonomy.parquet.snappy")
        tax_df.to_parquet(tax_path)

        os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        from .db_build_common import install_built_database
        return install_built_database(tmp_out, out_dir)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
