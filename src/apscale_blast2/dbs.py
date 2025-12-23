"""Database specification and validation utilities.

A database is considered valid when BLAST indices are present and the taxonomy
mapping table exists.
"""

from __future__ import annotations
import os
from dataclasses import dataclass
from .taxmap import find_taxmaps_paths

@dataclass
class DatabaseSpec:
    path: str  # database root folder or prefix

def is_blast_index_dir(d: str) -> bool:
    if not os.path.isdir(d): return False
    names = os.listdir(d)
    ok_idx = (any(n.endswith(".nin") or n.endswith(".nhr") for n in names) and any(n.endswith(".nsq") for n in names)) or any(n.endswith(".nal") for n in names)
    return ok_idx

def ensure_db_prefix(p: str) -> str:
    if os.path.isdir(p):
        if os.path.basename(p).lower()=="db" and is_blast_index_dir(p):
            for n in os.listdir(p):
                if n.endswith(".nal"): return os.path.join(p, n[:-4])
            for n in os.listdir(p):
                if n.endswith(".nin") or n.endswith(".nhr"): return os.path.join(p, n[:-4])
        cand = os.path.join(p, "db")
        if is_blast_index_dir(cand):
            for n in os.listdir(cand):
                if n.endswith(".nal"): return os.path.join(cand, n[:-4])
            for n in os.listdir(cand):
                if n.endswith(".nin") or n.endswith(".nhr"): return os.path.join(cand, n[:-4])
        if is_blast_index_dir(p):
            for n in os.listdir(p):
                if n.endswith(".nal") or n.endswith(".nin") or n.endswith(".nhr"):
                    return os.path.join(p, n.rsplit(".",1)[0])
    if os.path.exists(p+".nin") or os.path.exists(p+".nhr") or os.path.exists(p+".nal"):
        return p
    return p

def validate_database(db_path: str) -> None:
    pref = ensure_db_prefix(db_path)
    base_dir = os.path.dirname(pref) if os.path.isfile(pref) else (pref if os.path.isdir(pref) else os.path.dirname(pref))
    # indices
    idx_ok = False
    if os.path.isdir(base_dir):
        for n in os.listdir(base_dir):
            if n.endswith(".nal"):
                idx_ok = True; break
        if not idx_ok:
            have_nin = any(n.endswith(".nin") for n in os.listdir(base_dir))
            have_nhr = any(n.endswith(".nhr") for n in os.listdir(base_dir))
            have_nsq = any(n.endswith(".nsq") for n in os.listdir(base_dir))
            idx_ok = have_nsq and (have_nin or have_nhr)
    if not idx_ok:
        raise SystemExit(f"Invalid database: no BLAST indices were found in '{base_dir}'. Expected *.nal or *.nin/*.nhr + *.nsq")
    # taxonomy mapping file(s) near the database
    taxmaps = find_taxmaps_paths(pref)
    if not taxmaps:
        raise SystemExit(f"Invalid database: 'db_taxonomy.parquet.snappy' (or an equivalent CSV) was not found near '{pref}'.")
