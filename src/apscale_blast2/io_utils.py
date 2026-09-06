"""I/O helpers.

FASTA splitting and FASTA order preservation utilities.
"""

from __future__ import annotations
import os, pathlib, csv, gzip, re
from typing import List, Dict

FA_EXTS = (".fa", ".fasta", ".fna")
NUCLEOTIDES = re.compile(r"^[ACGTURYSWKMBDHVNacgturyswkmbdhvn]+$")


def fasta_stem(path: str) -> str:
    name = pathlib.Path(path).name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    return pathlib.Path(name).stem


def is_fasta(path: str) -> bool:
    return str(path).lower().endswith(FA_EXTS + tuple(x + ".gz" for x in FA_EXTS))


def iter_fasta(path: str):
    """Yield complete records and validate nucleotide FASTA structure."""
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    header, sequence = None, []
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    if not sequence:
                        raise ValueError(f"Empty FASTA sequence: {header.split()[0]}")
                    yield header, "".join(sequence)
                header, sequence = line[1:].strip(), []
                if not header:
                    raise ValueError(f"Empty FASTA header at line {number}: {path}")
            else:
                if header is None or not NUCLEOTIDES.fullmatch(line):
                    raise ValueError(f"Invalid nucleotide FASTA at line {number}: {path}")
                sequence.append(line)
        if header is None:
            raise ValueError(f"No FASTA records: {path}")
        if not sequence:
            raise ValueError(f"Empty FASTA sequence: {header.split()[0]}")
        yield header, "".join(sequence)


def read_fasta_order(path: str) -> List[str]:
    order, seen = [], set()
    for header, _ in iter_fasta(path):
        identifier = header.split()[0]
        if identifier in seen:
            raise ValueError(f"Duplicate FASTA ID {identifier!r}: {path}")
        seen.add(identifier)
        order.append(identifier)
    return order

def split_fasta(path: str, out_dir: str, subset_size: int) -> List[str]:
    if subset_size < 1:
        raise ValueError("subset_size must be positive")
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
        out_fh = open(out_path, "x", encoding="utf-8")
        return out_fh
    try:
        for header, sequence in iter_fasta(path):
            if out_fh is None or n_in_subset >= subset_size:
                out_fh = new_out()
            n_in_subset += 1
            out_fh.write(f">{header}\n{sequence}\n")
    finally:
        if out_fh:
            out_fh.close()
    return subset_paths

def read_db_map(csv_path: str) -> Dict[str,str]:
    mapping = {}
    with open(csv_path, newline="", encoding="utf-8-sig") as fh:
        for i, row in enumerate(csv.DictReader(fh)):
            fa = (row.get("fasta") or row.get("FASTA") or "").strip()
            db = (row.get("db") or row.get("DB") or "").strip()
            if not fa or not db or fa in mapping:
                raise ValueError(f"Row {i+2} in {csv_path} is invalid or duplicated (expected columns 'fasta' and 'db').")
            selected = pathlib.Path(db).expanduser()
            mapping[fa] = str(selected if selected.is_absolute() else pathlib.Path(csv_path).resolve().parent / selected)
    if not mapping:
        raise ValueError(f"Empty database mapping: {csv_path}")
    return mapping
