"""Per-database defaults for apscale_blast2.

Used by the interactive wizard to persist identity thresholds per database.

The defaults file is stored inside the database folder (next to the BLAST
index files), so it moves with the database when shared.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Tuple
from .filtering import DEFAULT_THRESHOLDS, thresholds_to_dict


DEFAULTS_FILENAME = "apscale_blast2_defaults.json"


@dataclass(frozen=True)
class DbDefaults:
    """Defaults attached to a DB folder."""

    thresholds: str


def defaults_path(db_prefix: str) -> Path:
    """Return the defaults file path for a given BLAST DB prefix."""
    source = Path(db_prefix).expanduser().absolute()
    if source.name == DEFAULTS_FILENAME:
        return source
    # A future builder destination is a folder; an existing index identifies a prefix.
    is_prefix = any(source.parent.glob(source.name + ".n*")) if not source.is_dir() else False
    db_dir = source.parent if is_prefix else source
    if db_dir.name == "db" and any(db_dir.parent.glob("db_taxonomy.*")):
        db_dir = db_dir.parent
    return db_dir / DEFAULTS_FILENAME


def load_db_defaults(db_prefix: str) -> Optional[DbDefaults]:
    """Load per-DB defaults from JSON if present."""
    p = defaults_path(db_prefix)
    if not p.exists():
        return None
    data: Any
    with p.open("r", encoding="utf-8") as fh:
        data = json.load(fh)
    # Minimal schema: {"thresholds": "97,95,90,87,85"}
    if not isinstance(data, dict):
        raise ValueError(f"Invalid database defaults object: {p}")
    thresholds = str(data.get("thresholds", "")).strip()
    if not thresholds:
        return None
    thresholds_to_dict(thresholds)
    return DbDefaults(thresholds=thresholds)


def save_db_defaults(db_prefix: str, thresholds: str) -> Path:
    """Save per-DB defaults (currently only thresholds)."""
    p = defaults_path(db_prefix)
    thresholds_to_dict(thresholds)
    payload = {
        "thresholds": thresholds,
    }
    with p.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return p


def get_thresholds_for_db(db_prefix: str, global_thresholds: str | None) -> Tuple[str, str, Path]:
    """Return thresholds for a DB, falling back to global defaults.

    Returns a tuple: (thresholds, source, defaults_path).
    - source is either "db" (loaded from apscale_blast2_defaults.json) or "global".
    - defaults_path is always returned to make saving straightforward.
    """
    p = defaults_path(db_prefix)
    if global_thresholds is not None:
        thresholds_to_dict(global_thresholds)
        return global_thresholds, "argument", p
    loaded = load_db_defaults(db_prefix)
    if loaded is None:
        return DEFAULT_THRESHOLDS, "global", p
    return loaded.thresholds, "db", p
