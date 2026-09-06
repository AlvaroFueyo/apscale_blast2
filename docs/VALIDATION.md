# Validation and Maintainer Review

Date: 2026-09-06. Baseline: commit `729a6fda8eb8cbbc0ef62db39c1d5ed9af5e24a9`
(version 1.2.1). Changes are unreleased. This review concerns the general-purpose
metabarcoding package, not a study-specific pipeline.

## Scope and Scientific Decisions

Reviewed the original 17 Python modules, CLI, builders, output handling, taxonomy
normalization, documentation and packaging. Implemented regression fixes in the
local repository without replacing the working bioinformatics installation.

Maintainer decisions incorporated:

- Compatible tied hits may contribute known ranks while empty ranks are ignored;
  genuine contradictions are flagged and truncate the consensus.
- cf./aff./hybrid or similarly uncertain references are not definite species;
  original labels and uncertainty are retained.
- DiatBarcode uses RCM by default; Kociolek is an explicit alternative.
- Conflicting DiatBarcode paths continue with conservative consensus and an audit,
  instead of choosing a parent arbitrarily.
- Native searches in this campaign use megablast only. The blastn option and legacy
  forced parameters remain, with their configuration tested without blastn searches.

## Corrected Failure Modes

| Area | Defect and implemented correction |
| --- | --- |
| CLI automation | Unconditional threshold input blocked redirected stdin; explicit DB modes no longer prompt. |
| Settings | DB defaults could override explicit flags or be saved in a shared parent; precedence and storage corrected. Decimal, ordered, finite percentages are validated. |
| Python API | Invalid workers/chunk sizes and arbitrary scratch deletion; runtime options validated and user-supplied scratch directories never deleted. |
| Consensus | Entirely blank intermediate rank stopped comparison, hiding lower conflicts; compatible rank consensus now continues across blanks. |
| Names | Qualifiers/hybrids could become confident binomials; preserved provenance and conservative rank limits. |
| Taxonomy lookup | Unrelated ancestor tables and conflicting duplicates could be selected silently; local unambiguous association and duplicate validation. |
| Reference IDs | Leading zeros were lost; PR2 taxonomy suffix punctuation could break mapping. IDs remain strings and newly built PR2/UNITE indices use stable keys. |
| Input validation | Duplicate/empty/bad query FASTA could fail only after expensive work; validation precedes search. |
| Search lifecycle | Work continued after failures and diagnostics were hidden; bounded scheduling, cancellation and native stderr reporting. |
| Windows | Native LMDB could not open Unicode DB paths; use ASCII short DB paths only for native tools. |
| Output safety | Previous tables or intermediate files could be destroyed; explicit overwrite, owned staging, locking and rollback on ordinary publication errors. |
| Large outputs | All hits accumulated in memory and Excel row limits failed late; chunked processing, typed Parquet and multiple Excel sheets. |
| Archive installation | ZIP path handling could escape the intended filename; reject unsafe members. Validate native precompiled indices before replacing an installed database. |
| Builders | Existing empty targets could gain an extra nested folder; staged installation publishes the intended root. |
| MIDORI2 | Species qualifiers were truncated; preserve complete taxon tokens after removing the taxid. |
| SILVA | Variable-depth lineages were assigned positional ranks; require an explicit rank map or normalized table. |
| trnL | Official headerless CRUX and comma CSV were not read correctly; support those formats and named Excel columns. |
| DiatBarcode | Old single-sheet reader failed on v16 and discarded missing-species rows; support both layouts, both trees, gap counts, conflicts and missing ancestors. |
| Documentation | Incorrect flags/extensions, overly broad format promises and defaults claims; corrected README and new current usage guide. |
| Packaging | Deprecated license table; SPDX license metadata and compatible setuptools minimum. Version left for maintainer selection. |

## Automated Tests

Local validation: Windows x64, Python 3.13.14, BLAST+ 2.17.0+, pandas 3.0.5,
openpyxl 3.1.5 and pyarrow 25.0.1. The regression/native suite has 68 tests.
Tests cover rank conflicts, blank ranks, identity boundaries, permutation invariance,
legacy dominance, uncertainty, leading-zero IDs, cache invalidation, malformed FASTA,
ZIP traversal, preserving prior databases, failed writes, output collisions,
cancellation, gzip queries, Excel text cells/sheet splitting and empty Parquet chunks.

The opt-in native tests build a gzipped PR2 fixture, install a real precompiled bundle,
and invoke the CLI with stdin closed and cp1252 output in paths containing spaces
and non-ASCII characters. They verify explicit 99% species thresholds, known and
no-hit sequences, retained intermediates and refusal to overwrite without consent.
A corrupt but structurally plausible index cannot replace the previous database.

Final local result: all 68 tests pass, including the two native tests. Ruff checks
for syntax/undefined-name errors and unused imports pass. `git diff --check` passes
with the repository's normal line-ending configuration. Source distribution and
wheel build successfully; all 20 wheel modules and the CLI parser were imported
from an isolated installation, without replacing the production environment.

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
```

To include native tests, set `BLAST_BIN` to the directory containing the BLAST+
executables. Without it, two native tests are explicitly skipped. CI is configured
for Windows/Linux with Python 3.10 and 3.13; those remote jobs have not been executed
as part of this local review. This is not a claim of tested Linux/macOS compatibility.

## Official-Reference Integration

Downloaded reference material from official publishers; MIDORI2 used the local
GB272 file supplied by the maintainer. Full downloads remain under ignored `dbs/`.
Release URLs and checksums are recorded in [reference_sources.json](validation/reference_sources.json).
Each integration fixture uses the first 100 complete records, in source order.
These are deterministic format tests, not representative ecological samples.

| Reference | Version and tested format | Source |
| --- | --- | --- |
| MIDORI2 | GB272 srRNA, BLAST FASTA, local file | [MIDORI2](https://www.reference-midori.info/) |
| PR2 | 5.1.1, UTAX FASTA | [Official release](https://github.com/pr2database/pr2database/releases/tag/v5.1.1) |
| UNITE | 19.02.2025, general fungi FASTA | [Publisher DOI](https://doi.org/10.15156/BIO/3301229) |
| SILVA | 138.2 SSURef NR99 plus tax_slv rank map | [Official archive](https://ftp.arb-silva.de/release_138_2/Exports/) |
| trnL | CRUX FASTA and headerless taxonomy | [CALeDNA reference links](https://ucedna.com/reference-databases-for-metabarcoding) |
| DiatBarcode | 16.3 XLSX, RCM and Kociolek | [Official dataset](https://doi.org/10.15454/TOMBYZ) |

SILVA 138.2 is a pinned test release, not a claim that it is the latest; release 144
has changed taxonomy conventions. Support for newer dialects requires additional
fixtures. See the publisher's [taxonomy documentation](https://www.arb-silva.de/documentation/silva-taxonomy).

Built seven sample databases (six recipes, two DiatBarcode trees). Ran 100 self
queries for each using default masking and a separate unmasked megablast run:
1,400 query executions in total. Validated query-row preservation, reference mapping
for every retained hit and near-full perfect matches for canonical sufficiently long
self queries. Full results are in [official_results.json](validation/official_results.json).

All seven build/search checks pass. Important observed results are not hidden:

- MIDORI2: one of 100 self queries has no retained hit both masked and unmasked. It
  contains 131 N bases in 989 nucleotides. Noncanonical references are not required
  by the fixture assertions to have a perfect full-length alignment.
- trnL: seven short queries (44-55 nt) have no hit with default masking; all obtain
  a hit in the unmasked comparison. Application defaults were not changed.
- Kociolek: some sample references remain unresolved because the supplied tree lacks
  their ancestry. Missing nodes are explicitly exported, not borrowed from RCM.
- Several tests reach the default 30-target ceiling; consensus is conditional on
  reported candidates and cannot establish that no additional tied taxa exist.
- Additional PR2, SILVA and DiatBarcode references contain disagreements between
  their genus field and species binomial. Those species are cleared and their source
  inconsistency is explicitly flagged, rather than published as a coherent lineage.

The entire DiatBarcode workbook was also parsed: 10,507 sequence records, 71 gap
characters removed from the search-copy sequences. RCM has zero contradictory or
missing-ancestor records in this check. Kociolek has 21 conflict-affected records
and 185 with missing ancestors. These counts concern sequence records, not the number
of unique conflicting tree nodes. Synthetic tests independently verify the conflict
algorithm; the first-100-record fixture alone does not cover every official conflict.

Reproduction, after placing the documented downloads under `dbs/official`:

```bash
python tools/prepare_reference_samples.py --count 100
python tools/test_official_references.py --blast-bin /path/to/blast/bin
```

No manual download is currently needed. A biological accuracy evaluation would still
benefit from independently identified mock-community sequences and an expected
taxonomic table. Self-search is not a sensitivity/specificity benchmark.

## Postprocessing Performance

Same synthetic 30-hit/query assignment workload, no native search. Measurements are
single local runs and include output writing; they do not represent whole-pipeline
speedups on large reference databases.

| Implementation | Queries / hits | Chunk size | Seconds | Peak RSS MiB |
| --- | --- | --- | --- | --- |
| Original 1.2.1 review | 2,000 / 60,000 | 2,000 | 64.48 | approximately 191 |
| Revised, one chunk | 2,000 / 60,000 | 2,000 | 3.99 | 168.38 |
| Revised, bounded chunks | 2,000 / 60,000 | 100 | 4.52 | 111.77 |
| Revised, larger query count | 10,000 / 300,000 | 100 | 20.92 | 116.11 |

New exports include extra audit/provenance fields. Original and revised measurements
also differ in fixed validation/instrumentation overhead; do not treat the ratios as
a statistically replicated benchmark. Taxonomy-loading memory is not modeled by
these small synthetic reference tables.

```bash
python -m pip install psutil
python tools/benchmark_postprocessing.py --queries 2000 --subset-size 100
python tools/benchmark_postprocessing.py --queries 10000 --subset-size 100
```

## Release Checklist and Limits

- Select a new release number in both `pyproject.toml` and `src/apscale_blast2/__init__.py`.
  Locally built 1.2.1-labelled artifacts are verification artifacts, not a new release.
- Replace the existing `authors = [{name = "You"}]` placeholder with the maintainer's
  preferred attribution before publishing.
- Run the added CI matrix remotely; native validation here is Windows/Python 3.13 only.
- Rebuild old PR2/UNITE/SILVA/DiatBarcode databases when adopting the corrected builders.
  Existing precompiled taxonomy is not silently migrated or repaired.
- Audit downstream readers: Excel may contain multiple sheets, extra columns are added,
  and stricter input/overwrite checks may intentionally reject older unsafe workflows.
- Preserve original reference releases and sidecars. Index size/mtime is recorded,
  not a cryptographic hash of every potentially multi-gigabyte BLAST volume.
- Full reference builds/searches, all supported OS/Python combinations, independent
  biological truth sets, and interrupted-process resume are not comprehensively tested.
- Memory remains proportional to taxonomy size and query-ID validation; there is no
  automatic resume and no power-failure-proof multi-file transaction.

No repository push, release publication, production replacement or study-specific
taxonomy rerun was performed by this review.
