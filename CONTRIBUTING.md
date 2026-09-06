# Contributing and Releasing

Keep the package general-purpose: do not hard-code a study's paths, markers,
thresholds, sample naming conventions or reference database.

## Repository Contents

Keep source code, tests, CI configuration, packaging metadata, maintained documents,
small public FASTA examples and sanitized validation summaries in Git. The versioned
v1.1.2 manuals are historical references, not the current manual.

Keep these local: environments, caches, downloaded databases, native executables,
BLAST indices, experimental inputs, assignments, intermediates, full CI logs,
credentials and signed download URLs. Use `dbs/` for reference/development material,
`data/` for study inputs and `results/` for analyses. These directories are ignored.
Do not copy private samples into the public `examples/` directory.

`.gitignore` deliberately does not hide every FASTA, JSON, CSV or XLSX. Small fixtures
may belong in the repository; extensions alone cannot distinguish them from private
data. Review new files before staging. Ignoring a file neither deletes it locally
nor removes an already tracked file or its history. Avoid `git add -f` for local data.

The source archive uses `MANIFEST.in` to include the maintained documentation,
versioned manuals, tests, tools and four public examples. The wheel contains the
Python package and its installation metadata, not databases or test results.
Check built archives too: Git ignore rules do not control package contents.

## Development Checks

Use an isolated Python environment, then:

```bash
python -m pip install -e . build ruff
python -m unittest discover -s tests -v
ruff check src tests tools --select E9,F63,F7,F82,F401
python -m build
```

The fast tests use synthetic fixtures and need neither BLAST nor reference downloads.
Two native tests are opt-in and otherwise reported as skipped. Set `BLAST_BIN` to
the directory containing `blastn`, `makeblastdb` and `blastdbcmd` (or their Windows
`.exe` versions) to include them. Native test searches use megablast only.

```bash
# Linux/macOS
export BLAST_BIN=/path/to/ncbi-blast/bin
python -m unittest discover -s tests -v
```

```powershell
# Windows PowerShell
$env:BLAST_BIN = 'C:\tools\ncbi-blast\bin'
python -m unittest discover -s tests -v
```

Optional official-reference and performance checks are documented in
[VALIDATION.md](docs/VALIDATION.md); they are not part of the lightweight CI jobs.
CI currently runs Windows and Linux on Python 3.10 and 3.13. A green job confirms
the checks that ran, not native BLAST or biological accuracy on every platform.

Use UTF-8 and LF for text. `.gitattributes` enforces LF except for Windows `.bat` and
`.cmd` scripts; `.editorconfig` communicates this to editors. Do not change global
Git settings to accommodate one checkout. On Windows, compare resolved paths or
file identity, not the spelling of long paths versus their 8.3 aliases.

## Before Publishing 2.0.0

1. Keep `pyproject.toml` and `src/apscale_blast2/__init__.py` at the same version.
   The chosen version is `2.0.0`; a Git commit title alone does not set this metadata.
2. Confirm the package author attribution, update the changelog/release date, and
   review README, [the current manual](docs/USAGE.md) and the validation report.
3. Run the checks above and confirm all remote CI jobs pass on the exact release
   commit. Do not describe untested Python or operating-system versions as verified.
4. Review the staged content and ignored paths:

   ```bash
   git status --short
   git diff --cached --stat
   git diff --cached --check
   git ls-files -ci --exclude-standard
   git check-ignore -v dbs/example/reference.fasta results/example/query_taxonomy.xlsx
   ```

   `git ls-files -ci --exclude-standard` should be empty. If it is not, review each
   tracked file before deliberately untracking it; do not delete local datasets.
   Review staged diffs for credentials, personal paths and unpublished sample IDs.
5. Build from a clean checkout of the intended commit. Inspect the wheel and source
   archive for unexpected data, caches, credentials and stale version labels. Do not
   publish older verification artifacts from `dbs/release_check/`.
6. Create the release/tag and publish only the verified artifacts. Document the
   migration requirements; previously built databases are not silently converted.
