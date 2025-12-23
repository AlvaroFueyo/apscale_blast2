"""I/O helpers.

FASTA splitting and FASTA order preservation utilities.
"""

from __future__ import annotations
import os, pathlib, csv
from typing import List, Dict

def read_fasta_order(path: str) -> List[str]:
    order = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                order.append(line[1:].strip().split()[0])
    return order

def split_fasta(path: str, out_dir: str, subset_size: int) -> List[str]:
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
    subset_paths: List[str] = []
    idx = 0
    n_in_subset = 0
    out_fh = None
    def new_out():
        nonlocal idx, out_fh, n_in_subset
        if out_fh: out_fh.close()
        idx += 1; n_in_subset = 0
        out_path = os.path.join(out_dir, f"subset_{idx}.fasta")
        subset_paths.append(out_path)
        out_fh = open(out_path, "w", encoding="utf-8")
        return out_fh
    out_fh = new_out()
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if n_in_subset >= subset_size:
                    out_fh = new_out()
                n_in_subset += 1
            out_fh.write(line)
    if out_fh: out_fh.close()
    return subset_paths

def read_db_map(csv_path: str) -> Dict[str,str]:
    mapping = {}
    with open(csv_path, newline="", encoding="utf-8") as fh:
        for i, row in enumerate(csv.DictReader(fh)):
            fa = (row.get("fasta") or row.get("FASTA") or "").strip()
            db = (row.get("db") or row.get("DB") or "").strip()
            if not fa or not db:
                raise SystemExit(f"Row {i+2} in {csv_path} is invalid (expected columns 'fasta' and 'db').")
            mapping[fa] = db
    return mapping
