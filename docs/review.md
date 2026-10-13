# Review the complete pipeline

`make review` verifies the frozen input tree, builds missing analysis products,
creates `results/review/index.html`, and enforces the audit gate. Open that HTML
file locally; its thumbnails and workflow graph work without a network.
Source links require the surrounding checkout. `make reproduce` also builds the
report as part of its default `all` target.

## What is checked

- Twelve input files and two directory trees must match the archived identity,
  including paths, sizes and membership. This preflight runs on every Make invocation.
- All 82 ledger rows are checked: values, confidence bounds, units, estimands,
  inferential units, source paths and archived panel references. Extra, missing
  or duplicate claims fail.
- Both sensitivity tables are checked cell by cell, including all 15 model settings.
- Complete regional and scanner prediction CSVs must match their archived
  SHA-256 fingerprints. This checks every prediction, fold and row identifier;
  it also deliberately detects serialization changes.
- All nine generated figures must provide both PDF and PNG, with valid hashes.
- All nine revised editorial figures (PDF/PNG/SVG), 41 panel-source mappings,
  historical recipes, retained sources and source-product checks are included.
  Tiny reference snapshots permit numerical comparison when source CSV/JSON
  bytes differ by platform but the values reproduce.

`reference/` comes from the archived c03b09b commit. The workflow never updates
these expectations. The numerical tolerance is rtol=1e-9, atol=1e-10; missing
values are accepted only where the reference has a missing value.

## Read failures

The report-writing job succeeds even when a comparison fails, preserving
`results/audit/reproduction_report.json`. A separate `audit` rule refuses to
create `results/audit/PASS`. `make review` renders the reviewer report before
running this gate. Input failures leave `results/provenance/input_validation.json`.
A green historical report is not evidence of a new run: use the fresh report,
logs, and `results/review/source_hashes.json` from the run being reviewed.

## Scientific code and automated checks

Serial and distributed scanner paths call `models.scanner.fit_outer_fold` for
nested selection, calibration and prediction. Regional and definition-sensitivity
models share `models.regional.regional_estimator`. Protein memberships remain
frozen; the readable fold manifest has one protein/count/fold tuple per line.
No new model, split design, estimand or baseline is introduced by this cleanup.

`make check` includes a synthetic nested scanner run through preparation, five
fold outputs and receipts, and reduction. It compares serial and distributed
predictions, checks protein separation, rejects missing folds, and tests that
changing held-out labels cannot affect selection, thresholds or predictions.
Another integration test executes the actual Snakemake gate and confirms that a
failed audit's JSON survives. GitHub Actions runs these checks on Linux and macOS.
A complete DAG dry-run uses empty input placeholders, without executing or bypassing
the production input checks. These are fast structural checks; the full frozen-data experiment remains a
separate CPU run because its complete inputs are distributed separately.

## Figures and captions

The report shows freshly generated analysis figures alongside the selected
11 October editorial exports in `assets/review_figures/`. Every revised panel
links to its source data or authored vectors. Historical plotting recipes are
inspection records, with machine-specific paths replaced by placeholders.
They are not advertised as portable rendering commands.

The revised exports are retained and verified, not freshly rendered by the
scientific DAG. Figure 1 A/E, S1 B/D, and S3 A–E retain authored vector analyses;
S1 E uses the retained Table S3. Their historical analysis states are not all
represented by the current simplified tables. Draft captions are included;
S3 still needs its original method-specific detail when integrated. Archived
ledger panel references are preserved rather than silently remapped to the
revised artwork. The older OGT identity results remain supporting results and
are not restored to the revised main Figure 3.

Review the input/audit commit first, the shared fitting/CI commit second, and
figure provenance/reporting third. Generated evidence can then be inspected
without obscuring the computational changes. DuckDB remains only at the two
CSV-to-Parquet boundaries; this cleanup adds no runtime dependency.
