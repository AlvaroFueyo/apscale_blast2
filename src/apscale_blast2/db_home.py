"""Database home directory utilities.

Databases are stored in a per-user data directory by default, with optional
overrides via the APSCALE_BLAST2_DB_HOME environment variable or the --db-home
CLI argument.
"""

from __future__ import annotations

import os
from pathlib import Path


ENV_DB_HOME = "APSCALE_BLAST2_DB_HOME"


def default_db_home() -> str:
    """Return a per-user, writable directory to store BLAST databases.

    We avoid installing databases inside the Python package directory because that
    location may be read-only (e.g., Program Files) or be wiped on upgrades.
    """
    # User override
    env = os.environ.get(ENV_DB_HOME)
    if env:
        return os.path.abspath(os.path.expanduser(env))

    # Platform defaults
    if os.name == "nt":
        base = os.environ.get("LOCALAPPDATA") or os.environ.get("APPDATA") or str(Path.home())
        return os.path.join(base, "apscale_blast2", "db")

    # XDG on Linux, Application Support on macOS
    xdg = os.environ.get("XDG_DATA_HOME")
    if xdg:
        return os.path.join(os.path.expanduser(xdg), "apscale_blast2", "db")
    return os.path.join(str(Path.home()), ".local", "share", "apscale_blast2", "db")


def get_db_home(db_home: str | None = None, ensure: bool = True) -> str:
    """Resolve the database home and optionally create it."""
    p = os.path.abspath(os.path.expanduser(db_home)) if db_home else default_db_home()
    if ensure:
        os.makedirs(p, exist_ok=True)
    return p


def db_folder_for_name(db_home: str, name: str) -> str:
    safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in name.strip())
    if not safe:
        raise ValueError("Empty database name.")
    return os.path.join(db_home, f"db_{safe}")
