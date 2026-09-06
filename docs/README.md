# Documentation

## Current Version: 2.0

- [User manual](USAGE.md): installation, CLI, database formats, assignment semantics,
  audit columns, migration and recovery.
- [Validation report](VALIDATION.md): tested environments, official reference fixtures,
  CI findings, benchmarks and remaining limits.
- [Changelog](../CHANGELOG.md): changes between versions.
- [Contributing and release checklist](../CONTRIBUTING.md): tests, repository contents
  and publication checks.

`USAGE.md` is the maintained manual. The README provides a shorter introduction;
there is no current v2.0 PDF. A PDF generated from an older version is not a substitute
for the current guide.

## Historical Version: 1.1.2

The [Markdown manual](apscale_blast2_manual_v1.1.2.md),
[Quarto source](apscale_blast2_manual_v1.1.2.qmd) and
[PDF](apscale_blast2_manual_v1.1.2.pdf) are retained for users reproducing old analyses.
They describe v1.1.2, not the v2.0 CLI, safety checks or assignment semantics.
The archived PDF has not been regenerated.

## Validation Data

Only small, sanitized provenance and result summaries belong in
[`validation/`](validation/). Full references, query-level results, downloaded CI
logs and developer environments stay outside Git, normally under ignored `dbs/`
or `results/`. Reference datasets retain their publishers' licences and citation
requirements; this project's MIT licence does not relicense those datasets.
