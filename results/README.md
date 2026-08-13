# Reproduced paper artifacts

Only reviewer-facing products are retained here:

- `analysis/`: transparent CSV and JSON outputs for each scientific analysis;
- `figures/`: regenerated PDF and PNG figures plus their checksum manifest;
- `tables/`: the manuscript result ledger in reviewable CSV and compressed Parquet;
- `provenance/`: hashes of immutable analysis inputs;
- `audit/`: the terminal numerical and figure-integrity report.

Restartable shards, logs, and benchmark files are generated under `work/`,
`logs/`, and `benchmarks/` during execution. Those directories are intentionally
ignored by version control and may be removed after a successful archived run.
