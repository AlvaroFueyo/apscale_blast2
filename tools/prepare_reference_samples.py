"""Extract small, deterministic fixtures from locally downloaded official releases.

Large source files and generated samples remain under ignored dbs/.
"""

import argparse
import gzip
import io
import itertools
import json
from pathlib import Path
import shutil
import sys
import tarfile
import zipfile

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from apscale_blast2.streaming import sha256_file


def records(handle):
    header, sequence = None, []
    for line in handle:
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(sequence)
            header, sequence = line[1:].strip(), []
        elif line.strip():
            sequence.append(line.strip())
    if header is not None:
        yield header, "".join(sequence)


def write_sample(name, handle, root, count):
    rows = list(itertools.islice(records(handle), count))
    with (root / f"{name}.fasta").open("w", encoding="utf-8") as dest:
        for header, sequence in rows:
            dest.write(f">{header}\n{sequence}\n")
    with (root / f"{name}_queries.fasta").open("w", encoding="utf-8") as dest:
        for index, (header, sequence) in enumerate(rows):
            dest.write(f">q{index:04d}\n{sequence}\n")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=100)
    args = parser.parse_args()
    if args.count < 1:
        parser.error("count must be positive")
    official, samples = ROOT / "dbs" / "official", ROOT / "dbs" / "test_samples"
    samples.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for name, source in [("midori2", ROOT / "dbs" / "MIDORI2_UNIQ_NUC_GB272_srRNA_BLAST.fasta"),
                         ("pr2", official / "pr2_5.1.1_UTAX.fasta.gz"),
                         ("silva", official / "silva_138.2_NR99.fasta.gz")]:
        opener = gzip.open if source.suffix == ".gz" else open
        with opener(source, "rt", encoding="utf-8") as handle:
            rows = write_sample(name, handle, samples, args.count)
        manifest[name] = {"source": source.name, "source_sha256": sha256_file(source), "records": len(rows)}
    shutil.copy2(official / "silva_138.2_ranks.txt.gz", samples / "silva_ranks.txt.gz")
    source = official / "unite_19.02.2025.tgz"
    with tarfile.open(source) as archive:
        member = next(m for m in archive.getmembers() if m.name.endswith("sh_general_release_dynamic_19.02.2025.fasta"))
        with archive.extractfile(member) as stream, io.TextIOWrapper(stream, encoding="utf-8") as handle:
            rows = write_sample("unite", handle, samples, args.count)
    manifest["unite"] = {"source": source.name, "member": member.name, "source_sha256": sha256_file(source), "records": len(rows)}
    source = official / "trnL_CRUX.zip"
    with zipfile.ZipFile(source) as archive:
        member = next(n for n in archive.namelist() if n.endswith("trnL_.fasta"))
        with archive.open(member) as stream, io.TextIOWrapper(stream, encoding="utf-8") as handle:
            rows = write_sample("trnl", handle, samples, args.count)
        ids = {header.split()[0] for header, _ in rows}
        tax_member = next(n for n in archive.namelist() if n.endswith("trnL_taxonomy.txt"))
        with archive.open(tax_member) as stream, io.TextIOWrapper(stream, encoding="utf-8") as handle, (samples / "trnl_taxonomy.txt").open("w", encoding="utf-8") as dest:
            for line in handle:
                if line.split("\t")[0] in ids:
                    dest.write(line)
    manifest["trnl"] = {"source": source.name, "member": member, "source_sha256": sha256_file(source), "records": len(rows)}
    source = official / "diatbarcode_16.3.xlsx"
    with pd.ExcelFile(source) as book, pd.ExcelWriter(samples / "diatbarcode.xlsx") as writer:
        sequences = pd.read_excel(book, sheet_name="sequences_info", dtype=str).fillna("").head(args.count)
        sequences.to_excel(writer, sheet_name="sequences_info", index=False)
        for sheet in ["taxo_RCM", "taxo_Kociolek"]:
            pd.read_excel(book, sheet_name=sheet, dtype=str).to_excel(writer, sheet_name=sheet, index=False)
    data = "".join(f">{identifier}\n{sequence.replace('-', '')}\n" for identifier, sequence in sequences[["Sequence ID", "Sequence"]].itertuples(index=False, name=None))
    write_sample("diatbarcode", io.StringIO(data), samples, args.count)
    manifest["diatbarcode"] = {"source": source.name, "source_sha256": sha256_file(source), "records": len(sequences)}
    for name in manifest:
        manifest[name]["sample_sha256"] = sha256_file(samples / f"{name}.fasta")
        manifest[name]["selection"] = "First N complete records, in source order; self-search tests, not a biological truth set"
    (samples / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
