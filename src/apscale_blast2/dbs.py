"""Database specification and validation utilities.

A database is considered valid when BLAST indices are present and the taxonomy
mapping table exists.
"""

from __future__ import annotations
from pathlib import Path
import re
from dataclasses import dataclass
from .taxmap import find_taxmaps_paths

@dataclass
class DatabaseSpec:
    path: str  # database root folder or prefix

def is_blast_index_dir(d: str) -> bool:
    return bool(_prefixes(Path(d)))


def _prefixes(directory: Path) -> list[Path]:
    if not directory.is_dir():
        return []
    aliases = sorted(directory.glob("*.nal"))
    candidates = [p.with_suffix("") for p in aliases]
    for index in sorted(directory.glob("*.nin")):
        prefix = index.with_suffix("")
        if all(Path(str(prefix) + ext).is_file() for ext in [".nin", ".nhr", ".nsq"]):
            # Volumes without an alias are selected by their shared parent prefix.
            prefix = Path(re.sub(r"\.\d+$", "", str(prefix)))
            if prefix not in candidates:
                candidates.append(prefix)
    return candidates

def ensure_db_prefix(p: str) -> str:
    path = Path(p).expanduser().absolute()
    if path.is_dir():
        candidates = _prefixes(path / "db") or _prefixes(path)
        if len(candidates) != 1:
            raise ValueError(f"Expected exactly one BLAST database in {path}; found {len(candidates)}")
        return str(candidates[0])
    if path not in _prefixes(path.parent):
        raise ValueError(f"Missing or incomplete BLAST database prefix: {path}")
    return str(path)

def validate_database(db_path: str) -> None:
    pref = ensure_db_prefix(db_path)
    # taxonomy mapping file(s) near the database
    taxmaps = find_taxmaps_paths(pref)
    if not taxmaps:
        raise ValueError(f"No taxonomy table associated with BLAST prefix: {pref}")
