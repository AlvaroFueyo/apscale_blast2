# Changelog

## 1.1.1 — Output polish, Parquet export, and runinfo fix

- Raw BLAST export: removed the redundant `Sequence ID` column.
- Raw BLAST export: `Accession` is now normalised (taxonomy suffixes such as `###...` are stripped; coordinate-like suffixes are reduced to `ACCESSION.VERSION` when possible).
- Added `--output-format {excel,parquet}` to write outputs as `.xlsx` (default) or compressed Parquet (`.parquet.snappy`).
- Fixed `runinfo` sidecar generation and ensured it is written to the `taxonomy/` output folder.
- Best-effort numeric dtype coercion for exported tables (keeps strings as strings, and metrics as int/float when possible).
- Fixed `Ambiguous taxa` formatting in `apscale2`: avoids duplicated genus in binomials (e.g., `Cyprinus carpio` instead of `Cyprinus Cyprinus carpio`) and filters rank placeholders like `genus`.

## 1.1.0 — New flag scheme (apscale2), legacy compatibility, and reproducibility

- Added `--flag-scheme` to select between:
  - `apscale2` (default): MRCA-based ambiguity handling with rank-aware trimming and Fl1/Fl2/... flags.
  - `apscale` (legacy): APSCALE-compatible behaviour including legacy defaults and legacy flag logic.
- Increased default `--max-target-seqs` to 30 in `apscale2` to better capture local database diversity.
- Added per-database thresholds support in the wizard via `apscale_blast2_defaults.json` (auto-created when user edits defaults).
- Added a `runinfo` sidecar (txt) capturing key run parameters (scheme, thresholds, task, max targets, DB metadata).
- Raw BLAST output now includes accession/sseqid plus mismatch and gap-open counts for easier manual verification.
- Bug fixes and UX polish in the wizard (e.g., database label shown with threshold prompts).

## 1.0.1 — Help output fix and improved readme

- Fixed `-h/--help` crash caused by unescaped percent signs in argparse help strings.
- Improved CLI help formatting with grouped options, examples, and a standard `-V/--version` flag.
- Improved Readme

## 1.0.0 — First public release

- Interactive wizard to select a FASTA folder and choose a database per FASTA.
- Local database store (per-user) with recipe-based builders:
  MIDORI2, UNITE, SILVA, PR2, trnL, and DiatBarcode.
- APSCALE-like hit selection (mode 1: Similarity → E-value) and F1–F4 flags.
- Query coverage strategy:
  - Hard filter in BLAST via -qcov_hsp_perc (default 50%)
  - Soft in-Python preference for >=90% (fallback to all hits if none pass)
- BLAST+ >= 2.17.0 requirement enforced at startup.
- In-memory taxonomy cache reused across multiple FASTA runs within one execution.

