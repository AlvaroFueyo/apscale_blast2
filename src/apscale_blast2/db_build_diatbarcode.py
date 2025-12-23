"""Database builder: DiatBarcode.

Builds a local BLAST database from the DiatBarcode Excel release.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from .db_home import db_folder_for_name


def build_diatbarcode_db(
    *,
    xlsx_path: str,
    db_home: str,
    name: str | None = None,
    makeblastdb_exe: str = "makeblastdb",
    keep_source: bool = True,
) -> str:
    """Build a DiatBarcode BLAST database from the official XLSX release.

    Replicates Till's script behaviour:
      - Reads sheet 'diatbarcode v12' (fallback: first sheet)
      - Writes a FASTA with 'Sequence ID' and 'Sequence'
      - Writes db_taxonomy.parquet.snappy with ranks extracted from the table
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

        # Read XLSX
        try:
            diat_df = pd.read_excel(p, sheet_name="diatbarcode v12").fillna("")
        except Exception:
            diat_df = pd.read_excel(p).fillna("")

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

        # Taxonomy table
        cols = [
            "Species",
            "Genus",
            "Family (following Round, Crawford & Mann 1990)",
            "Order (following Round, Crawford & Mann 1990)",
            "Class (following Round, Crawford & Mann 1990)",
            "Phylum (following Algaebase 2018)",
            "Subkingdom (following Algaebase 2018)",
            "Sequence ID",
        ]
        missing = [c for c in cols if c not in diat_df.columns]
        if missing:
            raise ValueError("Faltan columnas esperadas en el XLSX de DiatBarcode: " + ", ".join(missing))

        recs = []
        for row in diat_df[cols].values.tolist():
            species, genus, family, order, clazz, phylum, superkingdom, acc = [str(x).strip() for x in row]
            if not acc or not species:
                continue
            recs.append([acc, superkingdom, phylum, clazz, order, family, genus, species])
        tax_df = pd.DataFrame(
            recs,
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

        if keep_source:
            shutil.copy2(p, os.path.join(tmp_out, "source.xlsx"))

        # BLAST indices
        prefix = os.path.join(tmp_db_dir, "db")
        cmd = [makeblastdb_exe, "-in", fasta_path, "-title", "db", "-dbtype", "nucl", "-out", prefix]
        subprocess.run(cmd, check=True)

        tax_path = os.path.join(tmp_out, "db_taxonomy.parquet.snappy")
        tax_df.to_parquet(tax_path)

        os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        shutil.move(tmp_out, out_dir)
        return out_dir
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
