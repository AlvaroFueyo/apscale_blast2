"""Non-interactive database recipes for reproducible build scripts."""

import argparse

from .blast_tools import get_tool_version
from .db_defaults import save_db_defaults
from .db_home import get_db_home, install_precompiled_db_zip
from .filtering import thresholds_to_dict


def main(argv=None):
    parser = argparse.ArgumentParser(prog="apscale_blast2 build", description=__doc__)
    parser.add_argument("--recipe", required=True, choices=["midori2", "unite", "pr2", "silva", "trnl", "diatbarcode", "precompiled"])
    parser.add_argument("--input", required=True)
    parser.add_argument("--name")
    parser.add_argument("--db-home")
    parser.add_argument("--taxonomy", help="Accession/rank table (required for trnL; alternative to SILVA --rank-map).")
    parser.add_argument("--rank-map", help="Official SILVA tax_slv file matching the FASTA release, optionally gzipped.")
    parser.add_argument("--classification", choices=["RCM", "Kociolek"], default="RCM")
    parser.add_argument("--thresholds", help="Optional five identity thresholds saved with this database.")
    parser.add_argument("--makeblastdb-exe", default="makeblastdb")
    parser.add_argument("--blastdbcmd-exe", default="blastdbcmd", help="Native validation of precompiled database bundles.")
    parser.add_argument("--no-keep-source", action="store_true")
    args = parser.parse_args(argv)
    thresholds_to_dict(args.thresholds)
    if args.recipe == "trnl" and not args.taxonomy:
        parser.error("trnL requires --taxonomy")
    if args.recipe == "silva" and bool(args.taxonomy) == bool(args.rank_map):
        parser.error("SILVA requires exactly one of --taxonomy or --rank-map")
    if args.rank_map and args.recipe != "silva":
        parser.error("--rank-map is only valid for SILVA")
    if args.taxonomy and args.recipe not in {"silva", "trnl"}:
        parser.error("--taxonomy is only valid for SILVA or trnL")
    if args.classification != "RCM" and args.recipe != "diatbarcode":
        parser.error("--classification is only valid for DiatBarcode")
    if args.recipe != "precompiled" and (get_tool_version(args.makeblastdb_exe) or (0,)) < (2, 17, 0):
        parser.error("makeblastdb >= 2.17.0 is required")
    home = get_db_home(args.db_home)
    common = dict(db_home=home, name=args.name, makeblastdb_exe=args.makeblastdb_exe, keep_source=not args.no_keep_source)
    if args.recipe == "precompiled":
        built = install_precompiled_db_zip(args.input, home, args.name, blastdbcmd_exe=args.blastdbcmd_exe)
    elif args.recipe == "diatbarcode":
        from .db_build_diatbarcode import build_diatbarcode_db
        built = build_diatbarcode_db(xlsx_path=args.input, classification=args.classification, **common)
    elif args.recipe == "trnl":
        from .db_build_trnl import build_trnl_db
        built = build_trnl_db(input_fasta=args.input, taxonomy_path=args.taxonomy, **common)
    elif args.recipe == "silva":
        from .db_build_silva import build_silva_db
        built = build_silva_db(input_path=args.input, taxonomy_path=args.taxonomy, rank_map_path=args.rank_map, **common)
    else:
        from .db_build_midori2 import build_midori2_db
        from .db_build_pr2 import build_pr2_db
        from .db_build_unite import build_unite_db
        builder = {"midori2": build_midori2_db, "pr2": build_pr2_db, "unite": build_unite_db}[args.recipe]
        built = builder(input_path=args.input, **common)
    if args.thresholds:
        save_db_defaults(built, args.thresholds)
    print(f"Database installed: {built}", flush=True)
    return 0
