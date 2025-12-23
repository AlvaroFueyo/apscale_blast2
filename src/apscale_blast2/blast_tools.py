"""BLAST+ tool discovery and version checks.

apscale_blast2 requires NCBI BLAST+ >= 2.17.0 to ensure consistent behaviour
and to simplify database building workflows.
"""

from __future__ import annotations

import re
import subprocess


_VER_RE = re.compile(r"(\d+)\.(\d+)\.(\d+)")


def _parse_version(text: str) -> tuple[int, int, int] | None:
    m = _VER_RE.search(text or "")
    if not m:
        return None
    return int(m.group(1)), int(m.group(2)), int(m.group(3))


def get_tool_version(exe: str) -> tuple[int, int, int] | None:
    """Return (major, minor, patch) for a BLAST+ executable, or None if unknown."""
    try:
        r = subprocess.run([exe, "-version"], capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return None
    txt = (r.stdout or "") + "\n" + (r.stderr or "")
    return _parse_version(txt)


def require_blast_217(blastn_exe: str = "blastn", makeblastdb_exe: str = "makeblastdb") -> None:
    """Hard requirement: BLAST+ 2.17.0+ for portability (compressed FASTA support)."""
    v_blastn = get_tool_version(blastn_exe)
    v_mkdb = get_tool_version(makeblastdb_exe)
    if v_blastn is None:
        raise SystemExit(
            f"Not found '{blastn_exe}' en PATH. Install NCBI BLAST+ (>= 2.17.0) and try again."
        )
    if v_mkdb is None:
        raise SystemExit(
            f"Not found '{makeblastdb_exe}' en PATH. Install NCBI BLAST+ (>= 2.17.0) and try again."
        )

    if v_blastn < (2, 17, 0) or v_mkdb < (2, 17, 0):
        raise SystemExit(
            f"BLAST+ 2.17.0+ is required. Detected: blastn={v_blastn}, makeblastdb={v_mkdb}. "
            "Please upgrade BLAST+ and try again."
        )
