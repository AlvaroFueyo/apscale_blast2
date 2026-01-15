# apscale_blast2

`apscale_blast2` is a local BLAST-based taxonomic assignment tool inspired by [**apscale_blast**](https://github.com/TillMacher/apscale_blast).

## Project overview

Typical use case: metabarcoding workflows where you want to run local BLAST against
curated reference databases and obtain both **raw BLAST hits** and **taxonomy-aware assignments**.

## Key features

| Feature | apscale_blast | apscale_blast2 |
|--------|---------------|----------------|
| BLAST execution | Local + remote (NCBI) | Local only |
| Processing mode | Single FASTA per run | Batch processing of multiple FASTA files |
| Database reuse within run | No | Yes (taxonomy cached in memory) |
| Database handling | External, precompiled databases | Integrated database build and install |
| Assignment ranking | Similarity-first (mode 1) | Same (mode 1 replicated) |
| Taxonomic flags | F1–F4 | F1-F5 New logic (see below) |
| Query coverage handling | No | BLAST-level hard filter + soft post-filter |
| BLAST version requirement | Flexible | BLAST+ ≥ 2.17 required |
| Database location| User-defined path required for each run | Stored in local user data directory and auto-discovered |
| Output format | .xslx | .xslx or .parquet.snappy |

## Requirements

- Python **>= 3.10**
- [**NCBI BLAST+ >= 2.17.0**](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) in the Path (hard requirement). The program checks this at startup.

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

The wizard asks once for the BLAST search mode:
- `megablast` (default; faster; good for similar amplicons/barcodes)
- `blastn` (more sensitive; slower)

## CLI usage (non-interactive)

Run the same workflow without prompts by specifying inputs:

```bash
apscale_blast2 --fastas /path/to/fastas --db-for-all /path/to/db_folder --threads 8
```

Or provide a CSV mapping FASTA basenames to database folders:

```bash
apscale_blast2 --fastas /path/to/fastas --db-map mapping.csv
```

## Outputs

For each input FASTA, the tool writes:

- `raw_blast/<sample>_raw_blast.xlsx` — merged BLAST outfmt6 hits
- `taxonomy/<sample>_taxonomy.xlsx` — taxonomy-aware assignments after filtering/flags

Temporary subset FASTA files are created under a run directory and removed by default.

## Performance notes

- The tool applies early BLAST-side filters for speed:
  - `-evalue` (defaults to `1e-3`)
  - `-qcov_hsp_perc` (defaults to `50`)
- Additional **soft** query-coverage preference is applied in Python:
  - if any hit has `qcov >= 75%`, only those hits are considered
  - otherwise, the tool falls back to all hits (to avoid losing assignments)

## Reproducibility metadata (sidecar)

Each run writes a small sidecar file (e.g. `runinfo.txt`) alongside the outputs, including the effective parameters used
(`flag_scheme`, thresholds, BLAST task, `max_target_seqs`, and database path/name when available). This keeps tables clean
while still making results reproducible.

## Ambiguity flags

apscale_blast2 supports **two flagging/assignment schemes**, selectable via `--flag-scheme`.

### `--flag-scheme apscale2` (default)

Designed for **curated local databases** (including de-duplicated references) where “dominant taxon” heuristics are not very informative.

After applying identity/coverage thresholds and de-duplicating hits by taxon, the tool checks how many *unique taxa* remain and trims the final assignment to the **MRCA** (most recent common ancestor) when needed:

- **No flag:** only one species remains.
- **Fl1 — Two or more species (trimming to MRCA):** multiple species remain but they collapse to a single genus → final assignment is trimmed to the genus MRCA and the candidate species are stored under `Ambiguous taxa`.
- **Fl2/Fl3/… — Two or more genera/families/... (trimming to MRCA):** when multiple genera (or higher ranks) remain, the assignment is trimmed to the MRCA rank and all remaining candidates are stored under `Ambiguous taxa`.

If some hits are missing ranks (e.g. genus/species is empty), those missing values are ignored **when other hits provide a resolved value**, so you do not get false ambiguity just because one record is incompletely annotated.

### `--flag-scheme apscale` (legacy)

Replicates the **APSCALE / APSCALE-GUI** decision tree (F1–F4) described by Macher et al. (2023). When using this mode, apscale_blast2 also restores key legacy parameters (e.g. `task=blastn`, `max-target-seqs=20`, and no query-coverage filtering).

Macher T-H, Schütz R, Yildiz A, Beermann AJ, Leese F (2023) ﻿Evaluating five primer pairs for environmental DNA metabarcoding of Central European fish species based on mock communities. Metabarcoding and Metagenomics 7: e103856. [https://doi.org/10.3897/mbmg.7.103856](https://doi.org/10.3897/mbmg.7.103856)

## Wizard database defaults

When running the interactive wizard, apscale_blast2 looks for a per-database defaults file:

- `<db_folder>/apscale_blast2_defaults.json`

If present, the wizard loads and displays these defaults (identity thresholds, task, `max_target_seqs`, etc.) and lets you
edit them. If the file is missing, global defaults are shown; if you change any values, the wizard automatically creates
`apscale_blast2_defaults.json` inside the database folder so you do not need to remember them next time.

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

Apscale-blast2 could intall pre-compiled databases. The [pre-compiled databases are available under the following server](https://seafile.rlp.net/d/474b9682a5cb4193a6ad/) and will be updated regularly.


## License

MIT (see `LICENSE`).