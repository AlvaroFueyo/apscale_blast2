"""Per-database defaults for apscale_blast2.

Used by the interactive wizard to persist identity thresholds per database.

The defaults file is stored inside the database folder (next to the BLAST
index files), so it moves with the database when shared.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Tuple


DEFAULTS_FILENAME = "apscale_blast2_defaults.json"


@dataclass(frozen=True)
class DbDefaults:
    """Defaults attached to a DB folder."""

    thresholds: str


def defaults_path(db_prefix: str) -> Path:
    """Return the defaults file path for a given BLAST DB prefix."""
    db_dir = Path(os.path.dirname(db_prefix))
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
    thresholds = str(data.get("thresholds", "")).strip()
    if not thresholds:
        return None
    return DbDefaults(thresholds=thresholds)


def save_db_defaults(db_prefix: str, thresholds: str) -> Path:
    """Save per-DB defaults (currently only thresholds)."""
    p = defaults_path(db_prefix)
    payload = {
        "thresholds": thresholds,
    }
    with p.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return p


def get_thresholds_for_db(db_prefix: str, global_thresholds: str) -> Tuple[str, str, Path]:
    """Return thresholds for a DB, falling back to global defaults.

    Returns a tuple: (thresholds, source, defaults_path).
    - source is either "db" (loaded from apscale_blast2_defaults.json) or "global".
    - defaults_path is always returned to make saving straightforward.
    """
    p = defaults_path(db_prefix)
    loaded = load_db_defaults(db_prefix)
    if loaded is None:
        return global_thresholds, "global", p
    return loaded.thresholds, "db", p
