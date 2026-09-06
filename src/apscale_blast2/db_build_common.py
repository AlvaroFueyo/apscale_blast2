"""Shared, validated source extraction and database installation."""

from __future__ import annotations

import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import tempfile
import zipfile

from .io_utils import is_fasta, iter_fasta


def safe_member(name: str) -> None:
    normalized = name.replace("\\", "/")
    path = PurePosixPath(normalized)
    if path.is_absolute() or ".." in path.parts or ":" in normalized:
        raise ValueError(f"Unsafe archive member: {name!r}")


def resolve_input_file(input_path: str, workdir: str) -> str:
    source = Path(input_path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if source.suffix.lower() != ".zip":
        if not is_fasta(str(source)):
            raise ValueError(f"Expected nucleotide FASTA, FASTA.GZ or ZIP: {source}")
        return str(source)
    with zipfile.ZipFile(source) as archive:
        members = [entry for entry in archive.infolist() if not entry.is_dir()]
        for entry in members:
            safe_member(entry.filename)
        candidates = [entry for entry in members if is_fasta(entry.filename)]
        if len(candidates) != 1:
            raise ValueError("Source ZIP must contain exactly one FASTA; extract and select a file explicitly")
        entry = candidates[0]
        destination = Path(workdir).resolve() / Path(entry.filename.replace("\\", "/")).name
        with archive.open(entry) as src, destination.open("xb") as dst:
            shutil.copyfileobj(src, dst)
    return str(destination)


def install_built_database(source: str, destination: str) -> str:
    """Publish a completed build without nesting into an existing directory."""
    from .dbs import validate_database
    from .taxmap import load_taxmap_as_dict

    source_path, target = Path(source).resolve(), Path(destination).resolve()
    validate_database(str(source_path))
    load_taxmap_as_dict(str(source_path))
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if not target.is_dir() or any(target.iterdir()):
            raise FileExistsError(target)
    staging = Path(tempfile.mkdtemp(prefix=".apscale_blast2_build_", dir=target.parent))
    try:
        shutil.copytree(source_path, staging / "ready")
        if target.exists():
            target.rmdir()  # Fails safely if another process populated the folder.
        os.replace(staging / "ready", target)
    finally:
        shutil.rmtree(staging)
    return str(target)


def run_makeblastdb(executable: str, source: str, prefix: str, taxonomy) -> None:
    from .blast_runner import _resolve_sequence_id
    from .streaming import native_db_prefix
    from .taxmap import _normalize_columns

    normalized = _normalize_columns(taxonomy).fillna("")
    mapping = {}
    for row in normalized.itertuples(index=False, name=None):
        if row[0] in mapping and mapping[row[0]] != row[1:]:
            raise ValueError(f"Conflicting taxonomy for reference: {row[0]}")
        mapping[row[0]] = row[1:]
    # Use an uncompressed, validated source on BLAST's temporary filesystem.
    fasta = Path(prefix).parent.parent / "reference.fasta"
    seen = set()
    with fasta.open("x", encoding="utf-8") as handle:
        for header, sequence in iter_fasta(source):
            identifier = header.split()[0]
            key, ranks = _resolve_sequence_id({"sseqid": identifier}, mapping)
            if ranks is None:
                raise ValueError(f"Reference missing from taxonomy: {identifier}")
            if key in seen:
                raise ValueError(f"Duplicate reference FASTA identifier: {identifier}")
            seen.add(key)
            handle.write(f">{key}\n{sequence}\n")
    result = subprocess.run([executable, "-in", native_db_prefix(str(fasta)), "-title", "db", "-dbtype", "nucl", "-out", native_db_prefix(prefix)], capture_output=True, text=True, errors="replace")
    if result.returncode:
        raise RuntimeError(f"makeblastdb failed ({result.returncode}): {result.stderr.strip()}")
    fasta.unlink()
