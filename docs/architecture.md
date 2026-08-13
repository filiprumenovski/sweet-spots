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
which avoids OpenMP reduction-order drift. Lossless NumPy bundles carry model
features between preparation and fit stages; CSV remains the human-review
boundary for final tables. Kinase percentiles are reduced in original global
pair order, so changing the number of scheduling shards does not change sums.
