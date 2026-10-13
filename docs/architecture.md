# Workflow architecture

The workflow is a directed acyclic graph with one immutable source boundary and
one publication-facing audit boundary.

```text
frozen inputs + hashes
        |
  consensus regions
        |
        +-- self-clustering maps (PTM) --------------------+
        +-- peptide-aware null maps (protein shard) -------+
        +-- regional maps (outer fold, protein null) ------+
        +-- yin-yang maps (kinase peptide pair) -----------+--> figures
        +-- scanner maps (nested outer fold) --------------+--> metric ledger
        +-- evolution / OGT / FG-NUP ----------------------+--> final audit
```

## Map contract

A map job owns records by their stable global index. It writes its data to a
temporary file in the destination directory, flushes the file, and commits it
with an atomic rename. Only then does it write a JSON receipt containing shard
identity, record count, byte count, and SHA-256 for every output.

Shard count is an execution parameter, not an estimand. Random draws use global
record indices, never shard-local indices. Model jobs are partitioned by the
predeclared protein-grouped outer folds.

## Reduce contract

A reducer checks that receipts cover exactly `0..K-1`, that every output still
matches its receipted byte count and digest, that schemas agree, and that record
keys are non-null and unique. Fold reducers additionally require every prepared
row exactly once. Global operations such as false-discovery correction are
performed only after this validation.

## Failure model

- A killed writer cannot publish a truncated declared output.
- A completed shard is reusable after another worker fails.
- A modified or partially transferred shard fails its digest check.
- A changed shard count invalidates receipts rather than mixing partitions.
- A missing fold or duplicated record fails before publication-facing output.
- Logs and benchmark TSVs are independent per expensive job.

## Numerical determinism

Histogram-boosting jobs use one numerical thread. Parallelism is across folds,
which avoids OpenMP reduction-order drift. Protein-grouped train/test memberships
are frozen explicitly, including nested scanner folds, so CPU-dependent sorting
of equal-sized groups cannot change the experiment. The splitter rejects changed
cohorts. See [portable reproduction](portable_reproduction.md).

Lossless NumPy bundles carry model
features between preparation and fit stages; CSV remains the human-review
boundary for final tables. Kinase percentiles are reduced in original global
pair order, so changing the number of scheduling shards does not change sums.

## Review boundary

The immutable reference lives in `reference/`, outside writable `results/`.
Input identity is enforced before execution. The audit compares the full ledger,
sensitivity grids, prediction fingerprints, generated figures and retained
editorial provenance. The report job writes durable JSON; a separate gate fails
without deleting it. `make review` creates the offline HTML report before that
gate runs. See [the review guide](review.md).

All scientific rules declare the shared scientific source set, locked environment
and configuration as inputs. This conservatively rebuilds scientific jobs after
a scientific code change. Reporting code is tracked separately so presentation
changes do not retrain models. Source hashes in the report identify the files
used for the reviewed run.
