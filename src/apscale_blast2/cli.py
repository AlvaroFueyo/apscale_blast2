"""Command-line interface and interactive wizard.

This module implements:
- Non-interactive CLI execution (scripting/pipelines)
- An interactive wizard to select FASTA inputs and databases per FASTA
- Database recipe builders (MIDORI2, UNITE, SILVA, PR2, trnL, DiatBarcode)
"""

from __future__ import annotations
import argparse, os, shutil, logging
from typing import List, Tuple, Dict

from .blast_runner import run, RunOptions
from .dbs import DatabaseSpec, validate_database
from .io_utils import read_db_map
from .db_home import get_db_home, db_folder_for_name
from .db_build_midori2 import build_midori2_db
from .db_build_trnl import build_trnl_db
from .db_build_unite import build_unite_db
from .db_build_silva import build_silva_db
from .db_build_pr2 import build_pr2_db
from .db_build_diatbarcode import build_diatbarcode_db
from .blast_tools import require_blast_217

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
    p = argparse.ArgumentParser(prog="apscale_blast2",
        description=(
            "Select a folder of FASTA files, choose a database per FASTA, and run local BLAST (APSCALE-like). "
            "Databases are stored by default in a per-user local store (APSCALE_BLAST2_DB_HOME). "
            "Requires NCBI BLAST+ >= 2.17.0."
        ))
    p.add_argument("--fastas", help="Folder containing .fa/.fasta/.fna files (non-recursive). If omitted, you will be prompted.")
    p.add_argument("--db-home", help="Folder used to store/discover local databases (overrides the default).")
    p.add_argument("--dbs", help="Additional folder containing databases (optional). It will be added to the discovered list.")
    p.add_argument("--db-for-all", help="Use this database folder/prefix for all FASTA files (non-interactive mode).")
    p.add_argument("--db-map", help="CSV with columns 'fasta,db' (basenames and DB folders/prefixes). Non-interactive mode.")
    p.add_argument("--threads", type=int, default=0, help="Threads (0=auto=max(1,cpu-2)). Range: 0..1024 (default: 0).")
    p.add_argument("--workers", type=int, default=1, help="BLAST worker processes. Range: 1..128 (default: 1).")
    p.add_argument("--subset-size", type=int, default=100, help="Sequences per subset chunk. Range: 1..10000 (default: 100).")
    p.add_argument("--task", choices=["megablast","blastn"], default="megablast", help="def: megablast")
    p.add_argument("--max-target-seqs", type=int, default=20, help="Max hits per query. Range: 1..1000 (default: 20).")
    p.add_argument("--no-masking", action="store_true", help="Disable DUST/soft masking (default: enabled).")
    p.add_argument("--min-qcov", type=float, default=50.0, help="Minimum coverage % (qcovs/qcovhsp). Range: 0..100 (default: 50.0).")
    p.add_argument("--max-evalue", type=float, default=1e-3, help="Maximum E-value. Range: 1e-300..1 (default: 1e-3).")
    p.add_argument("--min-pident", type=float, default=50.0, help="Minimum identity % (pident). Range: 0..100 (default: 50).")
    p.add_argument("--inline-perc-identity", dest="inline_perc_identity", action="store_true", default=True,
                   help="Pass -perc_identity to BLAST using --min-pident (default: enabled).")
    p.add_argument("--no-inline-perc-identity", dest="inline_perc_identity", action="store_false",
                   help="Disable -perc_identity filtering during BLAST execution.")
    p.add_argument("--keep-tsv", action="store_true", help="Keep intermediate TSV files (default: cleaned up automatically).")
    p.add_argument("--blastn-exe", default="blastn", help="blastn executable (default: blastn).")
    p.add_argument("--makeblastdb-exe", default="makeblastdb", help="makeblastdb executable (default: makeblastdb).")
    p.add_argument("--thresholds", default="97,95,90,87,85", help="Thresholds for species,genus,family,order,class (APSCALE defaults).")
    p.add_argument("--log-level", choices=["DEBUG","INFO","WARNING","ERROR"], default="INFO", help="Logging level (default: INFO).")
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
    print("  [S] Skip this FASTA (do not analyse)")
    for i, (lab, path) in enumerate(dbs, start=1):
        print(f"  [{i}] {lab}")
    if last_choice:
        print(f"Press ENTER to reuse [{last_choice}]")
    while True:
        s = input("Choose DB (number / S): ").strip()
        if s == "" and last_choice is not None:
            return str(last_choice)
        if s.lower() in {"s", "skip"}:
            return "SKIP"
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

    logging.basicConfig(level=getattr(logging, a.log_level), format="%(levelname)s: %(message)s")
    logger = logging.getLogger("apscale_blast2")

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
    db_home = get_db_home(a.db_home, ensure=True)
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
            raise SystemExit("El CSV de --db-map no contiene estos FASTA: " + ", ".join(missing))
        selections = [(fa, mp[os.path.basename(fa)]) for fa in fastas]
        print("Non-interactive mode: using the --db-map mapping.", flush=True)
    else:
        if not dbs:
            print(f"\nDetected {len(fastas)} FASTA(s). No databases installed yet in: {db_home}", flush=True)
        else:
            print(f"\nDetected {len(fastas)} FASTA(s) y {len(dbs)} base(s) candidata(s).", flush=True)

        # Search-mode selection (interactive mode only)
        a.task = prompt_task_choice(default=a.task)

        last: int | None = None
        for fa in fastas:
            sel = choose_db_for_fasta(fa, dbs, last)
            if sel == "SKIP":
                selections.append((fa, None))
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
            print(f"  - {os.path.basename(fa)}  →  (saltado)", flush=True)
        else:
            print(f"  - {os.path.basename(fa)}  →  {dbp}", flush=True)

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

            print(f"  - Construyendo {recipe.upper()}: {name}", flush=True)
            if recipe == "trnl":
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

    opts = RunOptions(
        threads=int(a.threads), workers=int(a.workers), subset_size=int(a.subset_size),
        task=a.task, max_target_seqs=int(a.max_target_seqs), masking=(not a.no_masking),
        min_qcov=float(a.min_qcov), max_evalue=float(a.max_evalue), min_pident=float(a.min_pident),
        blastn_exe=a.blastn_exe, keep_tsv=bool(a.keep_tsv), thresholds=a.thresholds,
        log_level=a.log_level, inline_perc_identity=bool(a.inline_perc_identity),
    )

    tmp_root = os.path.join(fasta_dir, "_apscale_blast2_tmp")
    if os.path.exists(tmp_root): shutil.rmtree(tmp_root, ignore_errors=True)
    os.makedirs(tmp_root, exist_ok=True)

    for fa, dbp in selections:
        if dbp is None:
            print(f"Omitido: {os.path.basename(fa)}", flush=True)
            continue
        out_dir = os.path.join(tmp_root, os.path.splitext(os.path.basename(fa))[0])
        os.makedirs(out_dir, exist_ok=True)
        dbspec = DatabaseSpec(path=dbp)
        res = run(fa, out_dir, dbspec, opts)
        print(f"Listo: {res['filtered_xlsx']}", flush=True)

    try:
        if os.path.isdir(tmp_root) and not os.listdir(tmp_root):
            os.rmdir(tmp_root)
    except Exception:
        pass

    print("\nDone. Check the 'taxonomy/' folder in the parent directory of your FASTA files.", flush=True)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
