# apscale_blast2

`apscale_blast2` is a local BLAST-based taxonomic assignment tool inspired by **apscale_blast**.
It preserves the APSCALE-style assignment logic (mode-1 hit selection and F1–F4 flags) while adding
an interactive wizard and built-in database management.

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
| Taxonomic flags | F1–F4 | F1–F4 (unchanged logic) |
| Query coverage handling | No | BLAST-level hard filter + soft post-filter |
| BLAST version requirement | Flexible | BLAST+ ≥ 2.17 required |
| Database location| User-defined path required for each run | Stored in local user data directory and auto-discovered |


## Requirements

- Python **>= 3.10**
- [**NCBI BLAST+ >= 2.17.0**](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) in the Path (hard requirement). The program checks this at startup.

You can find the latest blast+ executables and further information on the installation [here](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. BMC Bioinformatics, 10, 421. [https://doi.org/10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)

## Installation

From the repository root:

```bash
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\Activate.ps1
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
  - if any hit has `qcov >= 90%`, only those hits are considered
  - otherwise, the tool falls back to all hits (to avoid losing assignments)

## Ambiguity flags (F1–F4)

apscale_blast2 follows the same *APSCALE/APSCALE-GUI* taxonomic assignment decision tree described by Macher et al. (2023), where ambiguous species-level BLAST results are handled using four flags (F1–F4). :contentReference[oaicite:0]{index=0}

### When are flags applied?
Flags are only relevant when **multiple species-level candidates remain** after the initial trimming / filtering steps (e.g., similarity-based ranking and threshold trimming). The goal is to either (i) keep a species call when it is clearly supported, or (ii) **downgrade** the assignment to genus or higher ranks when the ambiguity cannot be resolved. :contentReference[oaicite:1]{index=1}

### Flag definitions

- **F1 — Dominant taxon**  
  If, after trimming, a *dominant* species-level taxon remains, it is selected as the final assignment (“F1 – Dominant taxon”). :contentReference[oaicite:3]{index=3}

- **F2 — Two species, one genus**  
  If exactly **two species** of the **same genus** remain, the assignment is stored as the genus plus both possible species epithets separated by a slash (e.g., *Leuciscus idus/leuciscus*). :contentReference[oaicite:4]{index=4}

- **F3 — Multiple species of one genus**  
  If **one genus** remains but **more than two species** are still possible, the assignment is saved at genus level (e.g., *Hucho sp.*), with the list of ambiguous species retained as metadata. :contentReference[oaicite:5]{index=5}

- **F4 — Multiple genera (Trimming to MRCA)**  
  If candidates belong to **more than one genus** and no dominant taxon is present, the assignment is trimmed to the **most recent common taxon** (MRCA; Most Recent Common Ancestor, the lowest shared rank across the remaining candidates). :contentReference[oaicite:6]{index=6}

Macher T-H, Schütz R, Yildiz A, Beermann AJ, Leese F (2023) ﻿Evaluating five primer pairs for environmental DNA metabarcoding of Central European fish species based on mock communities. Metabarcoding and Metagenomics 7: e103856. [https://doi.org/10.3897/mbmg.7.103856](https://doi.org/10.3897/mbmg.7.103856)


## License

MIT (see `LICENSE`).

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