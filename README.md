# apscale_blast2

**Version 2.0** | [User manual](docs/USAGE.md) | [Validation](docs/VALIDATION.md) |
[Changelog](CHANGELOG.md) | [Contributing](CONTRIBUTING.md)

`apscale_blast2` is a local BLAST-based taxonomic assignment tool inspired by [**apscale_blast**](https://github.com/TillMacher/apscale_blast).

## Project overview

Typical use case: metabarcoding workflows where you want to run local BLAST against
curated reference databases and obtain both **raw BLAST hits** and **taxonomy-aware assignments**.

This is a general metabarcoding tool: it does not depend on an APSCALE project layout,
a particular marker, or a particular study. See the [current usage and interpretation guide](docs/USAGE.md)
and the [validation report](docs/VALIDATION.md). The [documentation index](docs/README.md)
distinguishes current instructions from the historical v1.1.2 manuals.

## Key features

| Feature | apscale_blast | apscale_blast2 |
|--------|---------------|----------------|
| BLAST execution | Local + remote (NCBI) | Local only |
| Processing mode | Single FASTA per run | Batch processing of multiple FASTA files |
| Database reuse within run | No | Yes (taxonomy cached in memory) |
| Database handling | External, precompiled databases | Integrated database build and install |
| Assignment ranking | Similarity-first (mode 1) | Same (mode 1 replicated) |
| Taxonomic flags | F1–F4 | Fl1-Fl7 (default); F1-F4 (legacy) |
| Query coverage handling | No | BLAST-level hard filter + soft post-filter |
| BLAST version requirement | Flexible | BLAST+ ≥ 2.17 required |
| Database location| User-defined path required for each run | Stored in local user data directory and auto-discovered |
| Output format | .xlsx | .xlsx or .parquet.snappy |

## Requirements

- Python **>= 3.10**
- [**NCBI BLAST+ >= 2.17.0**](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html),
  including `blastn` and `makeblastdb`; `blastdbcmd` is also needed to validate
  precompiled bundles. Put these on `PATH` or pass their executable paths.

Windows, Linux and macOS use the same Python package but need BLAST binaries for
their OS/architecture. The declared Python minimum is not a claim that every newer
version has been tested. See the [platform and CI evidence](docs/VALIDATION.md#remote-ci-log-review)
for verified environments and outstanding checks; macOS remains unverified.

You can find the latest blast+ executables and further information on the installation [here](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. BMC Bioinformatics, 10, 421. [https://doi.org/10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)

## Installation

You can install directly from the repository:

```bash
pip install "git+https://github.com/AlvaroFueyo/apscale_blast2.git"
```

If you are developing or modifying the code, a local editable install is more convenient:

```bash
git clone https://github.com/AlvaroFueyo/apscale_blast2.git
cd apscale_blast2
python -m venv venv

# Linux/macOS:
source venv/bin/activate
# Windows PowerShell:
venv\Scripts\Activate.ps1

pip install -U pip
pip install -e .
```

## Quick start (wizard mode)

Run:

```bash
apscale_blast2
```

You will be prompted to:

1. Select a folder containing `.fa/.fasta/.fna` files (non-recursive).
2. For each FASTA, choose:
   - **Build and install a new database (recipe)**, or
   - Select an existing installed database, or
   - **Skip** the FASTA.

The wizard asks once for the **BLAST search mode**:
- Search mode: `megablast` (default; faster; good for similar amplicons/barcodes) or `blastn` (more sensitive; slower)

The wizard uses the `apscale2` flagging scheme by default. The legacy `apscale` scheme remains available from the CLI via `--flag-scheme apscale`.

## CLI usage (non-interactive)

Run the same workflow without prompts by specifying inputs:

```bash
apscale_blast2 --fastas /path/to/fastas --db-for-all /path/to/db_folder --threads 8
```

Or provide a CSV mapping FASTA basenames to database folders:

```bash
apscale_blast2 --fastas /path/to/fastas --db-map mapping.csv
```

`--fastas` also accepts a single nucleotide FASTA, including `.fa.gz`, `.fasta.gz`,
and `.fna.gz`. Relative database paths in the mapping CSV are resolved relative to
the CSV, not the current working directory. Explicit database selection never asks
for input; `--db-for-all` and `--db-map` are mutually exclusive.

```bash
apscale_blast2 --fastas reads.fasta.gz --db-for-all dbs/db_reference --out-dir results/run1 --thresholds 99,95,90,87,85 --output-format parquet
apscale_blast2 build --recipe pr2 --input pr2_UTAX.fasta.gz --db-home dbs --name pr2
apscale_blast2 build --help
```

## Outputs

For each input FASTA, the tool writes:

- `raw_blast/<sample>_raw_blast.xlsx` — merged BLAST outfmt6 hits
- `taxonomy/<sample>_taxonomy.xlsx` — taxonomy-aware assignments after filtering/flags

Temporary subset FASTA files are created under a run directory and removed by default.

By default, output folders are placed in the parent of the FASTA directory. Use
`--out-dir` to select a different root. Existing results are protected; replacing
them requires `--overwrite`. Outputs are staged and published only after successful
processing. Excel files split across sheets before the row limit; Parquet is
recommended for large analyses. The historical Python `run(..., out_dir, ...)`
argument no longer controls scratch-file deletion; set `RunOptions.output_dir`.

## Performance notes

- The tool applies early BLAST-side filters for speed:
  - `-evalue` (defaults to `1e-3`)
  - `-qcov_hsp_perc` (defaults to `50`)
- Additional **soft** query-coverage preference is applied in Python:
  - if any hit has `qcov >= 75%`, only those hits are considered
  - otherwise, the tool falls back to all hits (to avoid losing assignments)

## Reproducibility metadata (sidecar)

Each input writes `<sample>.runinfo.json` and a legacy text sidecar alongside the taxonomy
table. JSON records effective options, BLAST version, counts, query/taxonomy SHA-256,
source-module hashes and index file metadata. Preserve these files with the results.
Raw hits also retain original species names, uncertainty and reference-mapping status.

## Ambiguity flags

apscale_blast2 supports **two flagging/assignment schemes**, selectable via `--flag-scheme`.

### `--flag-scheme apscale2` (default)

Designed for **curated local databases** (including de-duplicated references). The old dominance criterion was removed because apscale_blast2 works with local curated databases, usually already taxonomically deduplicated, where redundancy-driven “dominant taxon” heuristics lose much of their meaning.

After hit selection and rank-specific identity trimming, the tool compares unique
taxonomic profiles using a **rank-wise compatible consensus**. It does not calculate
an MRCA from an external taxonomy or phylogenetic tree:

- **No flag:** surviving trimmed lineages are compatible, ignoring missing ranks.
- **Fl1 — Two species of one genus / More than two species of one genus:** when all surviving taxa belong to the same genus, the final assignment keeps the genus and reports either `Genus epithet1/epithet2` (exactly two species) or `Genus sp.` (more than two species). The surviving taxa are listed under `Ambiguous taxa`.
- **Fl2-Fl7, conflicts from genus to kingdom:** the first conflicting rank and all
  lower ranks are cleared. Candidate profiles are retained under `Ambiguous taxa`.

If some hits are missing ranks (e.g. genus/species is empty), those missing values are ignored **when other hits provide a resolved value**, so you do not get false ambiguity just because one record is incompletely annotated.

Use `assigned_rank` and `assignment_status` to interpret results: Fl1 composite
labels are not resolved species, and cf./aff./hybrid names remain conservative.
The manual details [uncertainty and audit columns](docs/USAGE.md#outputs-and-audit-columns).

### `--flag-scheme apscale` (legacy)

Retains the **APSCALE / APSCALE-GUI** decision branch (F1-F4) described by Macher et al.
(2023) and restores legacy search parameters (`task=blastn`, `max-target-seqs=20`,
and no query-coverage filtering). Input validation and conservative name cleaning
are shared with v2.0; this is not a bit-for-bit replay of historical results.

Macher T-H, Schütz R, Yildiz A, Beermann AJ, Leese F (2023) ﻿Evaluating five primer pairs for environmental DNA metabarcoding of Central European fish species based on mock communities. Metabarcoding and Metagenomics 7: e103856. [https://doi.org/10.3897/mbmg.7.103856](https://doi.org/10.3897/mbmg.7.103856)

## Wizard database defaults

When running the interactive wizard, apscale_blast2 looks for a per-database defaults file:

- `<db_folder>/apscale_blast2_defaults.json`

Only identity thresholds are stored in this file. The precedence is explicit
`--thresholds`, then database defaults, then `97,95,90,87,85`. Decimals are supported;
values must be finite percentages ordered species >= genus >= family >= order >= class.
The wizard can edit and save thresholds; non-interactive runs never prompt or save them.

## Raw BLAST output

The raw BLAST table includes additional columns useful for manual review, such as:

- subject identifier / accession (`sseqid` / `saccver` depending on the BLAST installation)
- `mismatch`
- `gapopen`

## Databases

#### Midori2

Leray, M., Knowlton, N., & Machida, R. J. (2022). MIDORI2: A collection of quality controlled, preformatted, and regularly updated reference databases for taxonomic assignment of eukaryotic mitochondrial sequences. Environmental DNA, 4(4), 894–907. [https://doi.org/10.1002/edn3.303](https://doi.org/10.1002/edn3.303)

#### Unite

Nilsson, R. H., Larsson, K.-H., Taylor, A. F. S., Bengtsson-Palme, J., Jeppesen, T. S., Schigel, D., Kennedy, P., Picard, K., Glöckner, F. O., Tedersoo, L., Saar, I., Kõljalg, U., & Abarenkov, K. (2019). The UNITE database for molecular identification of fungi: Handling dark taxa and parallel taxonomic classifications. Nucleic Acids Research, 47(D1), Article D1. [https://doi.org/10.1093/nar/gky1022](https://doi.org/10.1093/nar/gky1022)

When using the all eukaryote database, please cite it as follows:

Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2024): UNITE general FASTA release for eukaryotes 2. Version 04.04.2024. UNITE Community. [https://doi.org/10.15156/BIO/2959335](https://doi.org/10.15156/BIO/2959335)

Includes global and 3% distance singletons.

When using the fungi database, please cite it as follows:

Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2024): UNITE general FASTA release for Fungi 2. Version 04.04.2024. UNITE Community. [https://doi.org/10.15156/BIO/2959333](https://doi.org/10.15156/BIO/2959333)

Includes global and 3% distance singletons.

#### SILVA

Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., Peplies, J., & Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools. Nucleic Acids Research, 41(D1), D590–D596. [https://doi.org/10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)

#### pr2

Guillou, L., Bachar, D., Audic, S., Bass, D., Berney, C., Bittner, L., Boutte, C., Burgaud, G., de Vargas, C., Decelle, J., del Campo, J., Dolan, J. R., Dunthorn, M., Edvardsen, B., Holzmann, M., Kooistra, W. H. C. F., Lara, E., Le Bescot, N., Logares, R., … Christen, R. (2013). The Protist Ribosomal Reference database (PR2): A catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy. Nucleic Acids Research, 41(Database issue), D597–D604. [https://doi.org/10.1093/nar/gks1160](https://doi.org/10.1093/nar/gks1160)

#### diat.barcode

Rimet, F., Gusev, E., Kahlert, M., Kelly, M. G., Kulikovskiy, M., Maltsev, Y., Mann, D. G., Pfannkuchen, M., Trobajo, R., Vasselon, V., Zimmermann, J., & Bouchez, A. (2019). Diat.barcode, an open-access curated barcode library for diatoms. Scientific Reports, 9(1), Article 1. [https://doi.org/10.1038/s41598-019-51500-6](https://doi.org/10.1038/s41598-019-51500-6)

#### CRUX
[trnl database](https://ucedna.com/reference-databases-for-metabarcoding)

Curd, E. E., Gold, Z., Kandlikar, G. S., Gomer, J., Ogden, M., O’Connell, T., Pipes, L., Schweizer, T. M., Rabichow, L., Lin, M., Shi, B., Barber, P. H., Kraft, N., Wayne, R., & Meyer, R. S. (2019). Anacapa Toolkit: An environmental DNA toolkit for processing multilocus metabarcode datasets. Methods in Ecology and Evolution, 10(9), 1469–1475. [https://doi.org/10.1111/2041-210X.13214](https://doi.org/10.1111/2041-210X.13214)

#### Precompiled databases

Precompiled bundles can be installed with `apscale_blast2 build --recipe precompiled`.
The [external database collection](https://seafile.rlp.net/d/474b9682a5cb4193a6ad/)
is separate from this code repository. Select the intended reference release and
verify its taxonomy format; installation does not silently repair old mappings.

## Upgrading and Contributing

Before replacing an older installation, read [the v2.0 migration notes](docs/USAGE.md#upgrading-from-1x).
Keep reference databases, experimental data and run outputs outside the public source
tree or in ignored local directories. Tests, small examples and sanitized validation
summaries remain public to make the software reviewable and reproducible. See
[CONTRIBUTING.md](CONTRIBUTING.md) for repository and release checks.


## License

MIT (see `LICENSE`).
