"""Command-line interface and interactive wizard.

This module implements:
- Non-interactive CLI execution (scripting/pipelines)
- An interactive wizard to select FASTA inputs and databases per FASTA
- Database recipe builders (MIDORI2, UNITE, SILVA, PR2, trnL, DiatBarcode)
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
import textwrap
from typing import Dict, List, Tuple

from .blast_runner import run, RunOptions
from .dbs import DatabaseSpec, validate_database
from .db_defaults import get_thresholds_for_db, save_db_defaults
from .io_utils import read_db_map
from .db_home import get_db_home, db_folder_for_name, install_precompiled_db_zip
from .db_build_midori2 import build_midori2_db
from .db_build_trnl import build_trnl_db
from .db_build_unite import build_unite_db
from .db_build_silva import build_silva_db
from .db_build_pr2 import build_pr2_db
from .db_build_diatbarcode import build_diatbarcode_db
from .blast_tools import require_blast_217
from . import __version__

FA_EXTS = (".fa", ".fasta", ".fna")

def discover_fastas(fasta_dir: str) -> List[str]:
    return [os.path.join(fasta_dir, n) for n in sorted(os.listdir(fasta_dir))
            if os.path.isfile(os.path.join(fasta_dir, n)) and n.lower().endswith(FA_EXTS) and not n.startswith("subset_")]

def discover_dbs(db_dir: str) -> List[Tuple[str, str]]:
    """Discover usable databases under a directory.

    We only list folders that validate as APSCALE BLAST2 databases (BLAST indices + taxonomy table).

    Notes:
      - We intentionally do NOT include the DB-home root itself (it is just a container).
      - We scan the directory itself and its immediate children.
    """
    if not db_dir or not os.path.isdir(db_dir):
        return []

    out: List[Tuple[str, str]] = []
    seen: set[str] = set()

    # Candidates: db_dir itself and immediate subdirectories
    cands = [os.path.abspath(db_dir)]
    try:
        for name in sorted(os.listdir(db_dir)):
            cand = os.path.join(db_dir, name)
            if os.path.isdir(cand):
                # Skip hidden/temp folders
                if name.startswith('.') or name.startswith('_'):
                    continue
                cands.append(os.path.abspath(cand))
    except Exception:
        pass

    for cand in cands:
        ab = os.path.abspath(cand)
        if ab in seen:
            continue
        seen.add(ab)
        try:
            validate_database(ab)
        except SystemExit:
            continue
        except Exception:
            continue

        bn = os.path.basename(ab)
        name = bn[3:] if bn.startswith('db_') else bn
        label = name
        out.append((label, ab))

    return out

def validate_range(name: str, val: float, lo: float, hi: float, integer=False):
    if integer and abs(val - int(val)) < 1e-9:
        val = int(val)
    if val < lo or val > hi:
        raise SystemExit(f"{name} out of range [{lo},{hi}]: {val}")
    return int(val) if integer else val

def build_parser():
    class _Fmt(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        """Help formatter combining defaults + multi-line descriptions."""

    description = (
        "apscale_blast2: local BLAST-based taxonomic assignment with APSCALE-like trimming and ambiguity flags.\n\n"
        "Interactive mode (wizard) is the default when you do not provide non-interactive database arguments.\n"
        "Databases are stored and auto-discovered from a per-user local store by default.\n"
        "NCBI BLAST+ >= 2.17.0 is required."
    )
    epilog = textwrap.dedent(
        """\
        Examples:
          # Interactive wizard (prompts for FASTA folder and DB selection)
          apscale_blast2

          # Wizard, but pre-select the FASTA folder
          apscale_blast2 --fastas path/to/fasta_dir

          # Non-interactive: use one database for all FASTA files
          apscale_blast2 --fastas path/to/fasta_dir --db-for-all path/to/db_folder

          # Non-interactive: per-FASTA mapping from a CSV (columns: fasta,db)
          apscale_blast2 --fastas path/to/fasta_dir --db-map mapping.csv
        """
    )

    p = argparse.ArgumentParser(
        prog="apscale_blast2",
        description=description,
        epilog=epilog,
        formatter_class=_Fmt,
    )

    p.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show program version and exit.",
    )

    # Note: argparse uses old-style string formatting for help messages.
    # Avoid raw percent signs in help strings (use "percent" or escape as "%%").

    # Input selection
    g_in = p.add_argument_group("Inputs")
    g_in.add_argument(
        "--fastas",
        help="Directory containing .fa/.fasta/.fna files (non-recursive). If omitted, you will be prompted in wizard mode.",
    )

    # Database selection/discovery
    g_db = p.add_argument_group("Databases")
    g_db.add_argument("--db-home", help="Directory used to store and auto-discover local databases (overrides the default).")
    g_db.add_argument("--dbs", help="Additional directory containing databases (optional). Added to the discovery list.")
    g_db.add_argument("--db-for-all", help="Use this database folder/prefix for all FASTA files (non-interactive mode).")
    g_db.add_argument(
        "--db-map",
        help="CSV mapping FASTA basenames to database folders/prefixes (columns: fasta,db). Non-interactive mode.",
    )
    g_db.add_argument(
        "--print-db-home",
        action="store_true",
        help="Print the resolved local database directory and exit.",
    )
    g_db.add_argument(
        "--open-db-home",
        action="store_true",
        help="Open the local database directory in your system file explorer and exit.",
    )

    # Execution/performance
    g_exec = p.add_argument_group("Execution")
    g_exec.add_argument(
        "--threads",
        type=int,
        default=0,
        help="BLAST threads (0=auto=max(1,cpu-2)).",
    )
    g_exec.add_argument(
        "--workers",
        type=int,
        default=1,
        help="BLAST worker processes (use 1 unless you know what you are doing).",
    )
    g_exec.add_argument("--subset-size", type=int, default=100, help="Sequences per subset chunk.")

    # BLAST settings
    g_blast = p.add_argument_group("BLAST")
    g_blast.add_argument(
        "--task",
        choices=["megablast", "blastn"],
        default="megablast",
        help="BLAST task: megablast is faster for similar barcodes/amplicons; blastn is more sensitive.",
    )
    g_blast.add_argument("--max-target-seqs", type=int, default=30, help="Maximum hits reported per query.")

    g_flags = p.add_argument_group("flags")
    g_flags.add_argument(
        "--flag-scheme",
        choices=["apscale2", "apscale"],
        default="apscale2",
        help=(
            "Which flagging/assignment scheme to use. "
            "'apscale2' is the new MRCA-based scheme (default). "
            "'apscale' replicates the original apscale_blast behaviour (legacy)."
        ),
    )
    g_blast.add_argument("--no-masking", action="store_true", help="Disable DUST/soft masking (default: enabled).")
    g_blast.add_argument("--blastn-exe", default="blastn", help="Path/name of the blastn executable.")
    g_blast.add_argument("--makeblastdb-exe", default="makeblastdb", help="Path/name of the makeblastdb executable.")

    # Filtering / trimming
    g_flt = p.add_argument_group("Filtering")
    g_flt.add_argument(
        "--min-qcov",
        type=float,
        default=50.0,
        help="Hard minimum query coverage (percent) applied at BLAST level (qcovs/qcovhsp).",
    )
    g_flt.add_argument(
        "--max-evalue",
        type=float,
        default=1e-3,
        help="Maximum E-value (passed to BLAST and used in post-processing).",
    )
    g_flt.add_argument(
        "--min-pident",
        type=float,
        default=50.0,
        help="Minimum percent identity (pident).",
    )
    g_flt.add_argument(
        "--inline-perc-identity",
        dest="inline_perc_identity",
        action="store_true",
        default=True,
        help="Pass -perc_identity to BLAST using --min-pident.",
    )
    g_flt.add_argument(
        "--no-inline-perc-identity",
        dest="inline_perc_identity",
        action="store_false",
        help="Do not pass -perc_identity to BLAST (apply identity filter only in post-processing).",
    )
    g_flt.add_argument(
        "--thresholds",
        default="97,95,90,87,85",
        help="Comma-separated thresholds for species, genus, family, order, class (APSCALE defaults).",
    )

    # Output/debug
    g_out = p.add_argument_group("Output")
    g_out.add_argument("--keep-tsv", action="store_true", help="Keep intermediate TSV files (default: cleaned up).")
    g_out.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging verbosity.",
    )

    return p

def prompt_dir(prompt_text: str) -> str:
    while True:
        p = input(prompt_text).strip().strip('"').strip("'")
        if p and os.path.isdir(p):
            return p
        print("Invalid path, please try again.")

def choose_db_for_fasta(fa: str, dbs: List[Tuple[str, str]], last_choice: int | None) -> str:
    print(f"\nFASTA: {os.path.basename(fa)}")
    print("  [0] Build and install a new database (recipe)")
    print("  [I] Install a precompiled database (.zip)")
    print("  [S] Skip this FASTA (do not analyze)")
    for i, (lab, path) in enumerate(dbs, start=1):
        print(f"  [{i}] {lab}")
    if last_choice:
        print(f"Press ENTER to reuse [{last_choice}]")
    while True:
        s = input("Choose DB (number / S / I): ").strip()
        if s == "" and last_choice is not None:
            return str(last_choice)
        if s.lower() in {"s", "skip"}:
            return "SKIP"
        if s.lower() in {"i", "install"}:
            return "INSTALL"
        if s.isdigit():
            k = int(s)
            if k == 0:
                return "0"
            if 1 <= k <= len(dbs):
                return str(k)
        print("Invalid input.")


def prompt_task_choice(default: str = "megablast") -> str:
    # Two options: megablast (fast) or blastn (more sensitive).
    print("\nBLAST search mode:")
    print("  [1] megablast  (faster; recommended for similar amplicons/barcodes)")
    print("  [2] blastn     (more sensitive; usually slower)")
    d = "1" if default == "megablast" else "2"
    while True:
        s = input(f"Choose mode (ENTER={d}): ").strip()
        if s == "":
            s = d
        if s in {"1","2"}:
            return "megablast" if s == "1" else "blastn"
        print("Invalid input.")


def prompt_builder_choice() -> str:
    print("\n— Build and install a new database —")
    opts = [
        ("midori2", "MIDORI2 (FASTA/FASTA.GZ/ZIP)") ,
        ("trnl", "trnL (FASTA + taxonomy TXT/CSV/TSV/XLSX)") ,
        ("unite", "UNITE (FASTA/FASTA.GZ/ZIP)") ,
        ("silva", "SILVA (FASTA; optional taxonomy table)") ,
        ("pr2", "PR2 (FASTA/FASTA.GZ/ZIP)") ,
        ("diatbarcode", "DiatBarcode (XLSX release)") ,
        ("precompiled_zip", "Precompiled BLAST database bundle (.zip; indices + taxonomy)") ,
    ]
    for i,(_,label) in enumerate(opts, start=1):
        print(f"  [{i}] {label}")
    while True:
        s = input("Choose recipe (number): ").strip()
        if s.isdigit() and 1 <= int(s) <= len(opts):
            return opts[int(s)-1][0]
        print("Invalid input.")


def _auto_name(recipe: str, src_path: str) -> str:
    base = os.path.basename(src_path)
    for suf in [".fasta.gz",".fa.gz",".fna.gz",".fasta",".fa",".fna",".gz",".zip",".xlsx",".xls"]:
        if base.lower().endswith(suf):
            base = base[: -len(suf)]
            break
    if base.lower().startswith("db_"):
        base = base[3:]
    base = base.replace(" ", "_")
    return f"{recipe}_{base}" if not base.lower().startswith(recipe) else base

def _prompt_file_path(prompt_text: str) -> str:
    while True:
        p = input(prompt_text).strip().strip('"').strip("'")
        if p and os.path.exists(p):
            return p
        print("Invalid path (does not exist), please try again.")

def _prompt_optional_name(prompt_text: str) -> str | None:
    s = input(prompt_text).strip()
    return s if s else None

def main(argv=None):
    p = build_parser()
    a = p.parse_args(argv)

    # Legacy APSCALE-BLAST mode: replicate the *behaviour* and the *defaults*
    # of the original apscale_blast implementation.
    # - task: blastn (not megablast)
    # - no query coverage filters
    if getattr(a, "flag_scheme", "apscale2") == "apscale":
        a.task = "blastn"
        a.min_qcov = 0.0
        a.max_target_seqs = 20

    logging.basicConfig(level=getattr(logging, a.log_level), format="%(levelname)s: %(message)s")
    logger = logging.getLogger("apscale_blast2")

    # Resolve DB home early so utility flags can work without requiring BLAST.
    db_home = get_db_home(a.db_home, ensure=True)

    if getattr(a, "print_db_home", False):
        print(db_home)
        return 0

    if getattr(a, "open_db_home", False):
        try:
            if os.name == "nt":
                os.startfile(db_home)  # type: ignore[attr-defined]
            elif sys.platform == "darwin":
                subprocess.run(["open", db_home], check=False)
            else:
                subprocess.run(["xdg-open", db_home], check=False)
        except Exception as e:
            raise SystemExit(f"Could not open database directory: {e}")
        return 0

    a.threads         = validate_range("--threads", float(a.threads), 0, 1024, integer=True)
    a.workers         = validate_range("--workers", float(a.workers), 1, 128, integer=True)
    a.subset_size     = validate_range("--subset-size", float(a.subset_size), 1, 10000, integer=True)
    a.max_target_seqs = validate_range("--max-target-seqs", float(a.max_target_seqs), 1, 1000, integer=True)
    a.min_qcov        = validate_range("--min-qcov", float(a.min_qcov), 0.0, 100.0)
    a.max_evalue      = validate_range("--max-evalue", float(a.max_evalue), 1e-300, 1.0)
    a.min_pident      = validate_range("--min-pident", float(a.min_pident), 0.0, 100.0)

    # Preflight BLAST tools (hard requirement)
    require_blast_217(a.blastn_exe, a.makeblastdb_exe)

    fasta_dir = a.fastas or prompt_dir("Folder with FASTA files: ")

    # Local DB store (no prompt): default per-user folder, overridable.
    search_dirs = [db_home]
    if a.dbs:
        if not os.path.isdir(a.dbs):
            raise SystemExit(f"--dbs is not a valid directory: {a.dbs}")
        search_dirs.append(a.dbs)

    fastas = discover_fastas(fasta_dir)
    if not fastas: raise SystemExit("No FASTA files were found.")
    dbs: List[Tuple[str,str]] = []
    for d in search_dirs:
        dbs.extend(discover_dbs(d))
    # dedup by absolute path
    seen=set(); dbs2=[]
    for lab,pp in dbs:
        ab=os.path.abspath(pp)
        if ab not in seen:
            seen.add(ab); dbs2.append((lab,pp))
    dbs = dbs2

    selections: List[Tuple[str, str | None]] = []
    # out_dir -> recipe + builder parameters (compiled up-front before running BLAST)
    build_reqs: Dict[str, Dict[str, str]] = {}
    if a.db_for_all:
        selections = [(fa, a.db_for_all) for fa in fastas]
        print("Non-interactive mode: using the same database for all FASTA files.", flush=True)
    elif a.db_map:
        from .io_utils import read_db_map
        mp = read_db_map(a.db_map)
        missing = [os.path.basename(fa) for fa in fastas if os.path.basename(fa) not in mp]
        if missing:
            raise SystemExit("The --db-map CSV does not contain these FASTA files: " + ", ".join(missing))
        selections = [(fa, mp[os.path.basename(fa)]) for fa in fastas]
        print("Non-interactive mode: using the --db-map mapping.", flush=True)
    else:
        if not dbs:
            print(f"\nDetected {len(fastas)} FASTA(s). No databases installed yet in: {db_home}", flush=True)
        else:
            print(f"\nDetected {len(fastas)} FASTA(s) and {len(dbs)} candidate database(s).", flush=True)

        # Search-mode selection (interactive mode only)
        a.task = prompt_task_choice(default=a.task)

        last: int | None = None
        for fa in fastas:
            sel = choose_db_for_fasta(fa, dbs, last)
            if sel == "SKIP":
                selections.append((fa, None))
                continue
            if sel == "INSTALL":
                src = _prompt_file_path("Path to the precompiled database ZIP (.zip): ")
                nm = _prompt_optional_name("Optional short name for the installed database (ENTER=auto): ")
                base = os.path.splitext(os.path.basename(src))[0]
                if base.lower().startswith("db_"):
                    base = base[3:]
                nm2 = nm or base
                safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in nm2.strip())
                if not safe.lower().startswith("db_"):
                    safe = "db_" + safe
                out_dir = os.path.join(db_home, safe)
                build_reqs[out_dir] = {"builder": "precompiled_zip", "input": src, "name": nm2}
                label = f"Precompiled:{os.path.basename(out_dir)}"

                # Make the installed database immediately available for the next FASTA files.
                dbs.append((label, out_dir))
                last = len(dbs)
                selections.append((fa, out_dir))
                continue

            if sel == "0":
                recipe = prompt_builder_choice()

                if recipe == "diatbarcode":
                    src = _prompt_file_path("Path to the DiatBarcode XLSX release: ")
                    nm = _prompt_optional_name("Short name for the database (ENTER=auto): ")
                    nm2 = nm or _auto_name(recipe, src)
                    out_dir = db_folder_for_name(db_home, nm2)
                    build_reqs[out_dir] = {"builder": recipe, "input": src, "name": nm2}
                    label = f"DiatBarcode:{os.path.basename(out_dir)}"
                elif recipe == "trnl":
                    src = _prompt_file_path("Path to the trnL FASTA file (.fasta/.fasta.gz/.zip): ")
                    tax = _prompt_file_path("Path to the trnL taxonomy file (TXT/CSV/TSV/XLSX): ")
                    nm = _prompt_optional_name("Short name for the database (ENTER=auto): ")
                    nm2 = nm or _auto_name(recipe, src)
                    out_dir = db_folder_for_name(db_home, nm2)
                    build_reqs[out_dir] = {"builder": recipe, "input": src, "taxonomy": tax, "name": nm2}
                    label = f"trnL:{os.path.basename(out_dir)}"
                elif recipe == "silva":
                    src = _prompt_file_path("Path to the SILVA FASTA file (.fasta/.fasta.gz/.zip): ")
                    tax = input("Optional path to a taxonomy table (ENTER=auto-parse from headers): ").strip().strip('"').strip("'")
                    if tax and not os.path.exists(tax):
                        print("Invalid path (does not exist). It will be ignored and headers will be parsed instead.")
                        tax = ""
                    nm = _prompt_optional_name("Short name for the database (ENTER=auto): ")
                    nm2 = nm or _auto_name(recipe, src)
                    out_dir = db_folder_for_name(db_home, nm2)
                    build_reqs[out_dir] = {"builder": recipe, "input": src, "taxonomy": tax, "name": nm2}
                    label = f"SILVA:{os.path.basename(out_dir)}"
                elif recipe == "unite":
                    src = _prompt_file_path("Path to the UNITE FASTA file (.fasta/.fasta.gz/.zip): ")
                    nm = _prompt_optional_name("Short name for the database (ENTER=auto): ")
                    nm2 = nm or _auto_name(recipe, src)
                    out_dir = db_folder_for_name(db_home, nm2)
                    build_reqs[out_dir] = {"builder": recipe, "input": src, "name": nm2}
                    label = f"UNITE:{os.path.basename(out_dir)}"
                elif recipe == "pr2":
                    src = _prompt_file_path("Path to the PR2 FASTA file (.fasta/.fasta.gz/.zip): ")
                    nm = _prompt_optional_name("Short name for the database (ENTER=auto): ")
                    nm2 = nm or _auto_name(recipe, src)
                    out_dir = db_folder_for_name(db_home, nm2)
                    build_reqs[out_dir] = {"builder": recipe, "input": src, "name": nm2}
                    label = f"PR2:{os.path.basename(out_dir)}"
                elif recipe == "precompiled_zip":
                    src = _prompt_file_path("Path to the precompiled database ZIP (.zip): ")
                    nm = _prompt_optional_name("Optional short name for the installed database (ENTER=auto): ")
                    base = os.path.splitext(os.path.basename(src))[0]
                    if base.lower().startswith("db_"):
                        base = base[3:]
                    nm2 = nm or base
                    safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in nm2.strip())
                    if not safe.lower().startswith("db_"):
                        safe = "db_" + safe
                    out_dir = os.path.join(db_home, safe)
                    build_reqs[out_dir] = {"builder": "precompiled_zip", "input": src, "name": nm2}
                    label = f"Precompiled:{os.path.basename(out_dir)}"
                else:
                    # midori2
                    src = _prompt_file_path("Path to the MIDORI2 file (.fasta/.fasta.gz/.zip): ")
                    nm = _prompt_optional_name("Short name for the database (ENTER=auto): ")
                    nm2 = nm or _auto_name(recipe, src)
                    out_dir = db_folder_for_name(db_home, nm2)
                    build_reqs[out_dir] = {"builder": recipe, "input": src, "name": nm2}
                    label = f"MIDORI2:{os.path.basename(out_dir)}"

                # Add the newly built database immediately (for the next FASTA files)
                dbs.append((label, out_dir))
                last = len(dbs)
                selections.append((fa, out_dir))
            else:
                k = int(sel)
                last = k
                selections.append((fa, dbs[k-1][1]))

    print("\nDatabase selection summary:", flush=True)
    for fa, dbp in selections:
        if dbp is None:
            print(f"  - {os.path.basename(fa)}  →  (skipped)", flush=True)
        else:
            print(f"  - {os.path.basename(fa)}  →  {dbp}", flush=True)

    # Per-database defaults (wizard): identity thresholds
    # We do this once per unique database so repeated FASTA→DB assignments
    # do not re-prompt.
    thresholds_cache: Dict[str, str] = {}
    for _, dbp in selections:
        if not dbp or dbp in thresholds_cache:
            continue
        current, src, defaults_path = get_thresholds_for_db(dbp, a.thresholds)
        # Human-friendly label for prompts/logging (database folder name).
        db_label = os.path.basename(os.path.normpath(dbp))
        print(f"\nIdentity thresholds for selected database: {db_label}", flush=True)
        print(f"  - Source: {src}", flush=True)
        print("  - Format: Species,Genus,Family,Order,Class", flush=True)
        print(f"  - Current: {current}", flush=True)
        new_val = input("  Enter new thresholds (or press ENTER to keep): ").strip()
        if new_val:
            # Validate early (and keep prompting until valid)
            while True:
                parts = [p.strip() for p in new_val.split(",") if p.strip()]
                ok = len(parts) == 5
                if ok:
                    try:
                        nums = [float(p) for p in parts]
                        ok = all(0.0 <= x <= 100.0 for x in nums)
                    except Exception:
                        ok = False
                if ok:
                    current = ",".join(parts)
                    try:
                        save_db_defaults(defaults_path, current)
                        print(f"  Saved: {defaults_path}", flush=True)
                    except Exception as e:
                        print(f"  Warning: could not save defaults file ({e}). Using the provided thresholds for this run.", flush=True)
                    break
                print("  Invalid value. Please enter 5 comma-separated numbers between 0 and 100.", flush=True)
                new_val = input("  Thresholds (Species,Genus,Family,Order,Class): ").strip()
                if not new_val:
                    break
        thresholds_cache[dbp] = current

    print("Preparing job: creating temporary files and configuring options...", flush=True)

    # 1) Build new databases (if any)
    if build_reqs:
        print("\nBuilding new databases (this may take a while)...", flush=True)
        for out_dir, info in build_reqs.items():
            # if already exists and looks valid, skip
            try:
                if os.path.isdir(out_dir) and os.listdir(out_dir):
                    validate_database(out_dir)
                    print(f"  - Already exists and is valid: {out_dir}", flush=True)
                    continue
            except Exception:
                pass
            recipe = info.get("builder") or "midori2"
            name = info.get("name") or os.path.basename(out_dir).replace("db_", "")
            src = info.get("input", "")
            tax = info.get("taxonomy", "")

            if recipe == "precompiled_zip":
                print(f"  - Installing PRECOMPILED: {name}", flush=True)
            else:
                print(f"  - Building {recipe.upper()}: {name}", flush=True)
            if recipe == "precompiled_zip":
                built = install_precompiled_db_zip(zip_path=src, db_home=db_home, name=name)
            elif recipe == "trnl":
                built = build_trnl_db(input_fasta=src, taxonomy_path=tax, db_home=db_home, name=name, makeblastdb_exe=a.makeblastdb_exe)
            elif recipe == "unite":
                built = build_unite_db(input_path=src, db_home=db_home, name=name, makeblastdb_exe=a.makeblastdb_exe)
            elif recipe == "silva":
                built = build_silva_db(input_path=src, taxonomy_path=(tax or None), db_home=db_home, name=name, makeblastdb_exe=a.makeblastdb_exe)
            elif recipe == "pr2":
                built = build_pr2_db(input_path=src, db_home=db_home, name=name, makeblastdb_exe=a.makeblastdb_exe)
            elif recipe == "diatbarcode":
                built = build_diatbarcode_db(xlsx_path=src, db_home=db_home, name=name, makeblastdb_exe=a.makeblastdb_exe)
            else:
                built = build_midori2_db(input_path=src, db_home=db_home, name=name, makeblastdb_exe=a.makeblastdb_exe)
            print(f"    OK → {built}", flush=True)

    # 2) Validate all selected databases
    from .dbs import validate_database
    for _, dbp in selections:
        if dbp is None:
            continue
        validate_database(dbp)

    # Base options (thresholds may vary per DB in wizard mode)
    opts_base = dict(
        threads=int(a.threads), workers=int(a.workers), subset_size=int(a.subset_size),
        task=a.task, max_target_seqs=int(a.max_target_seqs), masking=(not a.no_masking),
        min_qcov=float(a.min_qcov), max_evalue=float(a.max_evalue), min_pident=float(a.min_pident),
        flag_scheme=getattr(a, "flag_scheme", "apscale2"),
        blastn_exe=a.blastn_exe, keep_tsv=bool(a.keep_tsv),
        log_level=a.log_level, inline_perc_identity=bool(a.inline_perc_identity),
    )

    tmp_root = os.path.join(fasta_dir, "_apscale_blast2_tmp")
    if os.path.exists(tmp_root): shutil.rmtree(tmp_root, ignore_errors=True)
    os.makedirs(tmp_root, exist_ok=True)

    for fa, dbp in selections:
        if dbp is None:
            print(f"Skipped: {os.path.basename(fa)}", flush=True)
            continue
        thresholds = thresholds_cache.get(dbp, a.thresholds)
        opts = RunOptions(**opts_base, thresholds=thresholds)
        out_dir = os.path.join(tmp_root, os.path.splitext(os.path.basename(fa))[0])
        os.makedirs(out_dir, exist_ok=True)
        dbspec = DatabaseSpec(path=dbp)
        res = run(fa, out_dir, dbspec, opts)
        print(f"Done: {res['filtered_xlsx']}", flush=True)

    try:
        if os.path.isdir(tmp_root) and not os.listdir(tmp_root):
            os.rmdir(tmp_root)
    except Exception:
        pass

    print("\nDone. Check the 'taxonomy/' folder in the parent directory of your FASTA files.", flush=True)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
