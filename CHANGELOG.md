# Changelog

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

## 1.0.1 — Help output fix and improved readme

- Fixed `-h/--help` crash caused by unescaped percent signs in argparse help strings.
- Improved CLI help formatting with grouped options, examples, and a standard `-V/--version` flag.
- Improved Readme
