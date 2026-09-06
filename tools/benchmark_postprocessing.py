"""Measure post-processing only, with 30 synthetic hits/query and no BLAST search."""

import argparse
import json
from pathlib import Path
import sys
import tempfile
import threading
import time
from unittest.mock import patch

import psutil

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from apscale_blast2 import blast_runner as br
from apscale_blast2.dbs import DatabaseSpec
from apscale_blast2.io_utils import read_fasta_order


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queries", type=int, default=2000)
    parser.add_argument("--subset-size", type=int, default=100)
    args = parser.parse_args()
    work = Path(tempfile.mkdtemp(prefix="benchmark_", dir=ROOT / "dbs"))
    query = work / "benchmark.fasta"
    query.write_text("".join(f">q{i}\nACGTACGTACGTACGT\n" for i in range(args.queries)), encoding="ascii")
    database = work / "database"
    database.mkdir()
    for ext in ["nin", "nhr", "nsq"]:
        (database / f"db.{ext}").write_text("synthetic fixture, not a BLAST index", encoding="ascii")
    (database / "db_taxonomy.csv").write_text("Accession,superkingdom,phylum,class,order,family,genus,species\n" + "".join(f"ref{i},Eukaryota,Phyluma,Classa,Ordera,Familya,Genusa,Genusa alpha\n" for i in range(30)), encoding="ascii")

    def fake_blast(exe, prefix, fasta, *unused, **kwargs):
        out = Path(fasta + ".tsv")
        with out.open("w", encoding="ascii") as handle:
            for qid in read_fasta_order(fasta):
                for i in range(30):
                    handle.write(f"{qid}\tref{i}\tref{i}\tref{i}\t100\t1e-50\t100\t100\t0\t0\n")
        return str(out)

    process, stop = psutil.Process(), threading.Event()
    baseline = process.memory_info().rss
    peak = [baseline]

    def sample_memory():
        while not stop.wait(0.01):
            peak[0] = max(peak[0], process.memory_info().rss)

    sampler = threading.Thread(target=sample_memory)
    sampler.start()
    start = time.perf_counter()
    try:
        with patch.object(br, "_run_single", side_effect=fake_blast), patch("apscale_blast2.blast_tools.get_tool_version", return_value=None):
            br.run(str(query), None, DatabaseSpec(str(database)), br.RunOptions(threads=1, subset_size=args.subset_size, output_format="parquet", output_dir=str(work / "outputs")))
    finally:
        stop.set()
        sampler.join()
    result = {"queries": args.queries, "hits": args.queries * 30, "subset_size": args.subset_size, "seconds": time.perf_counter() - start,
              "baseline_rss_mib": baseline / 1024**2, "peak_rss_mib": peak[0] / 1024**2,
              "note": "Synthetic post-processing workload. Includes writing outputs; excludes BLAST searches and large database loading."}
    (work / "benchmark.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"Evidence: {work / 'benchmark.json'}")


if __name__ == "__main__":
    main()
