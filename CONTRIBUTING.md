# Contributing

Contributions that improve correctness, portability, or reviewability are
welcome. Please open an issue before changing an estimand, input cohort, null
model, or figure design so the scientific consequence is explicit.

## Development environment

```bash
git clone https://github.com/filiprumenovski/sweet-spots.git
cd sweet-spots
make install
make check
```

Python 3.12 is required. Dependencies are resolved only from the committed
`uv.lock`; workflow jobs invoke `uv run --frozen` and must not install packages
at runtime.

## Review contract

- Keep analysis logic in `src/regional_code_paper/analysis/` and orchestration
  in `workflow/`.
- Declare every input and output in the Snakemake DAG.
- Preserve protein grouping for resampling and cross-validation.
- Derive stochastic seeds through the shared randomness utilities.
- Write final tabular products as CSV or JSON; use the DuckDB CLI for Parquet
  materialization.
- Add a focused test for numerical or execution-contract changes.
- Do not commit provider-restricted data, credentials, logs, caches, or
  restartable workflow intermediates.
- Do not edit retained results by hand. Rebuild them through the workflow and
  update the audit in the same change.

## Before submitting a change

```bash
make check
make dry-run
```

If the change affects outputs, also run `make reproduce` and confirm that
`results/audit/reproduction_report.json` reports a successful reproduction.
