"""Database home directory utilities.

Databases are stored in a per-user data directory by default, with optional
overrides via the APSCALE_BLAST2_DB_HOME environment variable or the --db-home
CLI argument.
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile
from .db_build_common import safe_member


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

def install_precompiled_db_zip(zip_path: str, db_home: str, name: str | None = None, overwrite: bool = False, blastdbcmd_exe: str = "blastdbcmd") -> str:
    """Install a precompiled BLAST database bundle (.zip) into the local db home.

    The archive is expected to contain BLAST indices (e.g., .nsq/.nin/.nhr or .nal)
    and a taxonomy mapping table (e.g., db_taxonomy.parquet* or db_taxonomy.csv/tsv).

    Parameters
    ----------
    zip_path : str
        Path to the precompiled database archive (.zip).
    db_home : str
        Destination database home directory.
    name : str | None
        Optional database folder name. If omitted, the top-level folder inside the ZIP
        (if unambiguous) or the ZIP stem is used.
    overwrite : bool
        If True, overwrite an existing target directory.

    Returns
    -------
    str
        Path to the installed database directory.
    """
    zp = os.path.abspath(os.path.expanduser(zip_path))
    if not os.path.exists(zp):
        raise FileNotFoundError(zp)
    if not zp.lower().endswith(".zip"):
        raise ValueError("Precompiled database installer expects a .zip file.")

    os.makedirs(db_home, exist_ok=True)

    with zipfile.ZipFile(zp, "r") as zf:
        members = [m for m in zf.namelist() if not m.endswith("/")]
        if not members:
            raise ValueError("The ZIP archive is empty.")

        # Heuristic: detect BLAST indices inside the archive.
        idx_exts = (".nsq", ".nin", ".nhr", ".nal", ".ndb", ".not", ".ntf", ".nto")
        has_idx = any(m.lower().endswith(idx_exts) for m in members)
        if not has_idx:
            raise ValueError("No BLAST index files were found inside the ZIP archive.")

        # Derive a reasonable target folder name.
        top_levels = {m.split("/")[0] for m in members if "/" in m}
        if name:
            folder = name
        elif len(top_levels) == 1:
            folder = sorted(top_levels)[0]
        else:
            folder = os.path.splitext(os.path.basename(zp))[0]

    safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in folder.strip())
    if not safe:
        raise ValueError("Empty database name derived from archive.")
    if not safe.lower().startswith("db_"):
        safe = "db_" + safe

    target = os.path.join(db_home, safe)
    if os.path.exists(target) and not overwrite:
        raise FileExistsError(f"Target database folder already exists: {target}")

    tmp = tempfile.mkdtemp(prefix=".apscale_blast2_install_", dir=db_home)
    try:
        with zipfile.ZipFile(zp, "r") as zf:
            for member in zf.infolist():
                safe_member(member.filename)
            zf.extractall(tmp)

        # If the ZIP has a single top-level directory, install that directory; otherwise
        # install the extracted tree as-is into the target.
        candidates = [os.path.join(tmp, d) for d in os.listdir(tmp)]
        root_dir = candidates[0] if len(candidates) == 1 and os.path.isdir(candidates[0]) else tmp

        from .dbs import validate_database, ensure_db_prefix
        from .taxmap import load_taxmap_as_dict
        from .streaming import native_db_prefix
        validate_database(root_dir)
        load_taxmap_as_dict(root_dir)
        check = subprocess.run([blastdbcmd_exe, "-db", native_db_prefix(ensure_db_prefix(root_dir)), "-info"], capture_output=True, text=True, errors="replace")
        if check.returncode:
            raise ValueError(f"Cannot open the precompiled BLAST database: {check.stderr.strip()}")
        backup = None
        if os.path.exists(target):
            backup = tempfile.mkdtemp(prefix=".apscale_blast2_previous_", dir=db_home)
            os.rmdir(backup)
            os.replace(target, backup)
        try:
            os.replace(root_dir, target)
        except BaseException:
            if backup:
                os.replace(backup, target)
            raise
        if backup:
            shutil.rmtree(backup)
    finally:
        # If we moved tmp itself, it no longer exists; ignore errors.
        shutil.rmtree(tmp, ignore_errors=True)

    return target
