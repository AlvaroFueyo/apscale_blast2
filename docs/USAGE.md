# apscale_blast2 2.0 User Manual

This is the maintained manual for the 2.0 release line, not the historical v1.1.2 PDF.
The package is intended for general metabarcoding projects and standalone FASTA files.
It assigns taxonomy to nucleotide queries; read trimming, denoising, clustering and
negative-control removal belong to the upstream processing workflow.

## Installation and Requirements

Use Python >= 3.10 in a virtual environment. For a source checkout:

```bash
python -m venv .venv
```

Activate with `source .venv/bin/activate` on Linux/macOS, or
`.\.venv\Scripts\Activate.ps1` in Windows PowerShell, then install:

```bash
python -m pip install .
apscale_blast2 --version
apscale_blast2 --help
apscale_blast2 build --help
```

Pip installs the Python dependencies (pandas, openpyxl, pyarrow and tqdm), but not
NCBI BLAST+. Install BLAST+ >= 2.17.0 separately for the OS/architecture in use.
The assignment CLI checks `blastn` and `makeblastdb`; precompiled bundle validation
also needs `blastdbcmd`. Each can be supplied by its `--*-exe` option instead of
changing `PATH`. Verify executables with `blastn -version`, `makeblastdb -version`
and `blastdbcmd -version`; these commands do not run a search.

The Python package is shared across Windows, Linux and macOS, but native dependencies
must match the platform. See [the validation report](VALIDATION.md#remote-ci-log-review)
for actual test coverage. macOS and newer Python versions are not verified by that
report. No global Python or PATH change is required by these installation steps.

## Reproducible CLI

```bash
apscale_blast2 --fastas queries.fasta.gz --db-for-all dbs/db_reference --out-dir results/run1 --threads 8 --thresholds 99,95,90,87,85 --output-format parquet
```

The executable is `apscale_blast2` (underscore); the distribution installed by pip
is `apscale-blast2`. `python -m apscale_blast2.cli` is equivalent. The default task
is megablast. `--task blastn` remains available, but no blastn-task searches were run
in the current validation campaign. Legacy `--flag-scheme apscale` still forces
blastn, 20 targets and zero hard/soft coverage thresholds.

Quote paths containing spaces. FASTA IDs are the first whitespace-delimited token;
they must be unique within each input. Input sequences must be nonempty, unaligned
IUPAC nucleotide strings (RNA U is supported). Empty files, duplicate IDs and invalid
bases are rejected before BLAST. Gzipped nucleotide FASTA is supported. Protein
FASTA, FASTQ and arbitrary multi-member source archives are not accepted as queries.

Use one FASTA or a non-recursive directory. Files sharing a stem, such as `a.fa`
and `a.fasta.gz`, cannot share an output root. Separate experiments with the same
input filename should use different `--out-dir` roots.

`--db-map mapping.csv` accepts `fasta,db` columns; filenames include their extensions.
Relative database paths are relative to the CSV. The explicit database options do
not prompt and do not create the default user database store.

## Building Databases

All recipes require the exact documented dialect, not any file bearing a database name.
Use `apscale_blast2 build --help` for flags. Builders validate reference IDs, nucleotide
content, taxonomy coverage and conflicting duplicate IDs before publishing a build.
An existing nonempty destination is never overwritten by a recipe.

| Recipe | Supported input | Important details |
| --- | --- | --- |
| `midori2` | MIDORI2 BLAST FASTA with `###root_` and seven ranks | Preserves cf./aff./hybrid qualifiers before conservative cleaning. |
| `pr2` | Official UTAX FASTA (`;tax=k:...,p:...`) | Uses explicit UTAX rank codes, not taxo-long position guessing. Stable IDs exclude the taxonomy suffix. |
| `unite` | General FASTA with pipe-separated fields and `k__`, `p__`, etc. | Choose the intended normal/developer release explicitly after extracting TGZ. Stable IDs keep the first four pipe fields. |
| `silva` | FASTA plus matching official `tax_slv` map, or a normalized accession/rank table | Never infers ranks from positions in the semicolon lineage. |
| `trnl` | FASTA plus CRUX headerless taxonomy, named CSV/TSV or XLSX | CRUX format is accession TAB seven semicolon-separated ranks. Named tables use accession plus seven ranks. |
| `diatbarcode` | Legacy v12 flat XLSX or v16 sequence/tree XLSX | RCM default; `--classification Kociolek` selects the other tree. |
| `precompiled` | ZIP containing BLAST nucleotide indices and an associated taxonomy table | Uses `blastdbcmd -info` before publication. |

Plain FASTA, gzipped FASTA and ZIP with exactly one FASTA are accepted by FASTA
recipes. Extract TAR/TGZ or ambiguous ZIPs yourself and select the correct file;
the tool will not choose between biological releases silently. Archive paths with
traversal, absolute paths or drive prefixes are rejected.

```bash
apscale_blast2 build --recipe silva --input silva.fasta.gz --rank-map tax_slv_ssu_release.txt.gz --db-home dbs --name silva
apscale_blast2 build --recipe trnl --input trnL.fasta --taxonomy trnL_taxonomy.txt --db-home dbs --name trnl
apscale_blast2 build --recipe diatbarcode --input diatbarcode.xlsx --classification Kociolek --db-home dbs --name diatoms_kociolek
apscale_blast2 build --recipe precompiled --input reference.zip --db-home dbs --name reference
```

`--makeblastdb-exe` and `--blastdbcmd-exe` accept executable paths. `--thresholds`
on the build command optionally saves per-database thresholds. Sources are retained
by default; `--no-keep-source` disables copying them into the built database.

For DiatBarcode, gaps (`-`) and sequence-cell whitespace are removed from the copy
used to build the nucleotide index. Counts are recorded in `build_audit.csv`; the
original XLSX is not modified. Conflicting tree paths are reconciled rank by rank,
stopping at the first contradiction. All lower ranks are cleared. Missing ancestors
are not invented: known ranks remain, and missing nodes are recorded separately.
Tree cycles remain fatal. The audit includes original names and alternative lineages.
The workbook may include several markers; selecting a classification does not select
a marker. Filter the source deliberately when a marker-specific database is needed.

Database layout is normally `db_NAME/db/db.*` and `db_NAME/db_taxonomy.parquet.snappy`.
Custom databases may use an unambiguous taxonomy CSV/TSV/Parquet next to their index.
`db_taxonomy.*` is preferred; multiple candidate tables at the same level are rejected.
The loader does not borrow unrelated tables from ancestor directories.

## Assignment Semantics

The retained pipeline is: hard identity/coverage/e-value filters, soft coverage
preference (if any hit qualifies), maximum identity, minimum e-value among ties,
rank-specific identity trimming, then consensus across remaining taxonomic profiles.
Thresholds are species, genus, family, order, class. They are configurable heuristics,
not universally calibrated biological confidence levels. A 99% identity threshold
does not establish that a marker distinguishes two species.

For default `apscale2`, repeated reference entries do not outvote other taxa. Missing
rank values are ignored where other selected hits supply compatible information.
An entirely empty intermediate rank does not stop comparison of lower ranks.
The first real conflict truncates the assignment; no inferred ancestry is inserted.
This is a rank-wise compatible consensus, not an external taxonomy-tree or phylogenetic
MRCA calculation. Disjoint incomplete annotations are not proof of shared ancestry.

Names containing cf., aff., nr., group/complex qualifiers, unresolved species labels,
multiple names or hybrid notation are not converted into definite species. Original
species strings and uncertainty are retained. Hybrids do not yield a definite genus.
Uncertainty in one hit does not erase compatible positive information in another.
If a reference's binomial contradicts its own genus column, the species is cleared
and the inconsistency is flagged; the original label is preserved for review.

`Fl1` labels such as `Testus alpha/beta` or `Testus sp.` are display labels for
ambiguity, not species IDs. Use `assigned_rank` when counting species assignments.
Conflicts at genus through kingdom produce `Fl2` through `Fl7`. The legacy `apscale`
scheme preserves its dominance decision branch but shares the safer input and name
cleaning rules; it is not a bit-for-bit replay of every historical release.

## Outputs and Audit Columns

`raw_blast/` contains hits after hard filters, not every possible hit in the database.
`taxonomy/` contains one row per input query, in input order, including queries with
no retained hit. Both outputs support Excel or compressed Parquet.

| Column | Meaning |
| --- | --- |
| `assignment_status` | `assigned`, `ambiguous`, `unresolved`, `no_match`, or `taxonomy_missing`. |
| `assigned_rank` | Deepest resolved rank; does not count Fl1 composite labels as species. |
| `hit_count` | Hits retained after hard filtering, before soft preference and best-hit selection. |
| `candidate_taxa` | Unique selected taxonomic profiles after identity trimming. |
| `hit_limit_reached` | Retained distinct subjects reached `max_target_seqs`; additional candidates may be missing. |
| `taxonomy_mapped` | Raw-hit reference ID was found in the associated taxonomy table. |
| `species_original`, `species_uncertainty` | Raw reference name and conservative uncertainty category. |
| `reference_uncertainty` | Uncertainty categories among selected hits. |
| `reference_taxonomy_conflict` | Audited source-tree conflicts among selected references. |
| `reference_missing_nodes` | Missing source-tree ancestors, separate from contradictory ranks. |

Completely unmapped nonempty BLAST results cause failure rather than publishing a
misleading table. Partial mapping failures are counted, warned about and flagged.
If best selected hits include an unmapped reference, the status is `taxonomy_missing`;
known taxa may still be shown but should not be treated as an unqualified assignment.

Excel splits into `part_1`, `part_2`, etc. at 1,048,575 data rows per sheet. Readers must
read all sheets. Cells over 32,767 characters cause an explicit error directing users
to Parquet; values are not silently truncated. IDs remain text, including leading
zeros and strings beginning with `=`. Large Excel workbooks may still be impractical.

## Safety, Resources and Recovery

Jobs use unique output-side staging directories. Existing tables require explicit
`--overwrite`. A failed search or write leaves previous published outputs untouched;
partial diagnostics remain in the named work directory. `--keep-tsv` preserves only
this run's intermediates and records their location in the JSON sidecar.

Publication is atomic per file, with rollback on ordinary exceptions. It is not a
filesystem-wide transaction under abrupt power loss. The JSON sidecar is published
last; do not consume outputs while the same input's `.apscale_blast2.lock` exists.
After a crash, verify that no process owns the lock before manually removing it.
There is no automatic chunk resume yet; rerun into a new output root for recovery.

Only a bounded number of BLAST workers and hit chunks are in flight. Total requested
threads are divided across workers. Taxonomy itself is still loaded in memory;
memory also scales with query-ID validation and the chosen chunk size. Very large
reference mappings can still require substantial RAM. Bounded hit processing is not
a promise of constant memory for every component.

Windows BLAST/LMDB can fail on Unicode database paths. The runner uses an ASCII short
path for the external executable while retaining the original path for taxonomy.
If short names are unavailable, move the database to an ASCII path. This does not
alter global PATH or Python installations.

Default masking can remove short/low-complexity queries even on self-search. Treat
`no_match` as a result of the configured search, not evidence that the taxon is absent.
`max_target_seqs=30` is not an exhaustive census of tied references. Inspect the
hit-limit indicator and compare higher limits where taxonomic decisions depend on it.

## Upgrading from 1.x

Keep the original software environment, reference releases, settings and results
when reproducibility of an earlier analysis matters. Install v2.0 in a separate
environment and use a new output root for comparison.

- Rebuild databases made with older PR2, UNITE, SILVA or DiatBarcode builders to
  adopt corrected IDs and taxonomy handling. Installing an old precompiled ZIP does
  not convert its contents to the corrected interpretation.
- Replace interactive automation workarounds with explicit `--fastas`, `--db-for-all`
  or `--db-map`, and `--out-dir`. CLI thresholds now override database defaults.
- Use `assigned_rank` and `assignment_status`, not a nonempty `Species` cell alone,
  to count resolved species. Conservative name cleaning and conflict handling can
  intentionally produce different assignments from earlier versions.
- Read every Excel sheet or use Parquet. Additional audit columns are intentional;
  downstream code must not depend on a fixed column count or one-sheet workbooks.
- Existing outputs now require `--overwrite`. Failed runs retain diagnostics; they
  do not support automatic resume. Preserve JSON sidecars alongside final tables.
- Python callers should set `RunOptions.output_dir`. The historical positional
  `run(..., out_dir, ...)` parameter no longer defines a directory that can be deleted.

Do not put study inputs, full reference releases or run outputs into public commits.
See [the contribution guide](../CONTRIBUTING.md) for local directory conventions and
release checks, and [the changelog](../CHANGELOG.md) for the full change summary.
