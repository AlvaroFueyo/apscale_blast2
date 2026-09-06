"""Build and self-search samples of official databases with real BLAST+.

Run prepare_reference_samples.py first. Generated results stay in ignored dbs/.
Self-search validates file compatibility and traceability, not biological accuracy.
"""

import argparse
from dataclasses import replace
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import tempfile
import time
import traceback

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from apscale_blast2.blast_runner import RunOptions, run
from apscale_blast2.db_build_diatbarcode import build_diatbarcode_db
from apscale_blast2.db_build_midori2 import build_midori2_db
from apscale_blast2.db_build_pr2 import build_pr2_db
from apscale_blast2.db_build_silva import build_silva_db
from apscale_blast2.db_build_trnl import build_trnl_db
from apscale_blast2.db_build_unite import build_unite_db
from apscale_blast2.dbs import DatabaseSpec
from apscale_blast2.io_utils import read_fasta_order, iter_fasta
from apscale_blast2.taxmap import load_taxmap_as_dict


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast-bin", type=Path, required=True)
    args = parser.parse_args()
    sys.stdout.reconfigure(encoding="utf-8", errors="backslashreplace")
    samples = ROOT / "dbs" / "test_samples"
    destination = Path(tempfile.mkdtemp(prefix="validation_", dir=ROOT / "dbs"))
    executable = lambda name: str(args.blast_bin.resolve() / (name + (".exe" if sys.platform == "win32" else "")))
    report = {"date": datetime.now(timezone.utc).isoformat(), "sample_manifest": json.loads((samples / "manifest.json").read_text()), "results": {}, "output_root": str(destination)}
    for name in ["midori2", "pr2", "unite", "silva", "trnl", "diatbarcode", "diatbarcode_Kociolek"]:
        started = time.monotonic()
        try:
            recipe = name.split("_")[0]
            common = dict(db_home=str(destination), name=name, makeblastdb_exe=executable("makeblastdb"), keep_source=False)
            if recipe == "diatbarcode":
                built = build_diatbarcode_db(xlsx_path=str(samples / "diatbarcode.xlsx"), classification="Kociolek" if name.endswith("Kociolek") else "RCM", **common)
            elif recipe == "silva":
                built = build_silva_db(input_path=str(samples / "silva.fasta"), rank_map_path=str(samples / "silva_ranks.txt.gz"), **common)
            elif recipe == "trnl":
                built = build_trnl_db(input_fasta=str(samples / "trnl.fasta"), taxonomy_path=str(samples / "trnl_taxonomy.txt"), **common)
            else:
                builder = {"midori2": build_midori2_db, "pr2": build_pr2_db, "unite": build_unite_db}[recipe]
                built = builder(input_path=str(samples / f"{recipe}.fasta"), **common)
            query = samples / f"{recipe}_queries.fasta"
            options = RunOptions(threads=2, workers=2, subset_size=25, blastn_exe=executable("blastn"), output_dir=str(destination / name), output_format="parquet", thresholds="99,95,90,87,85")
            result = run(str(query), None, DatabaseSpec(built), options)
            taxonomy, raw = pd.read_parquet(result["taxonomy_output"]), pd.read_parquet(result["raw_output"])
            expected = len(read_fasta_order(str(query)))
            assert len(taxonomy) == expected, "Lost query rows"
            assert raw["taxonomy_mapped"].all(), "Unmapped reference hits"
            defaults_statuses = taxonomy["assignment_status"].value_counts().to_dict()
            # DUST can hide short/low-complexity self queries. Test unmasked megablast
            # separately, without changing the application's default masking policy.
            unmasked = run(str(query), None, DatabaseSpec(built), replace(options, masking=False, output_dir=str(destination / f"{name}_unmasked")))
            plain = pd.read_parquet(unmasked["taxonomy_output"])
            plain_raw = pd.read_parquet(unmasked["raw_output"])
            assert plain_raw["taxonomy_mapped"].all(), "Unmapped unmasked reference hits"
            assert not plain["assignment_status"].eq("taxonomy_missing").any(), "Missing taxonomy for an unmasked self match"
            clean_ids = {header for header, sequence in iter_fasta(str(query)) if set(sequence.upper()) <= set("ACGTU") and len(sequence) >= 60}
            near_full = plain_raw[(plain_raw["Similarity"] == 100) & (plain_raw["query_coverage"] >= 95)]
            assert clean_ids <= set(near_full["unique ID"]), "A canonical self query has no near-full perfect match"
            mapped = load_taxmap_as_dict(built)
            report["results"][name] = {"passed": True, "queries": expected, "hits": len(raw), "mapped_references": len(mapped), "default_statuses": defaults_statuses, "unmasked_statuses": plain["assignment_status"].value_counts().to_dict(), "canonical_self_queries_checked": len(clean_ids), "elapsed_seconds": time.monotonic() - started}
        except Exception as error:
            report["results"][name] = {"passed": False, "error": str(error), "traceback": traceback.format_exc()}
        print(json.dumps({name: report["results"][name]}, ensure_ascii=True), flush=True)
        (destination / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Evidence: {destination / 'report.json'}", flush=True)
    return 0 if all(r["passed"] for r in report["results"].values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
