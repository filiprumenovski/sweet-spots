# Immutable scientific reference

Copied from archived commit `c03b09bcb170767f8319bf08df18e1bedbed979e`, not from a new run.
The ledger covers all claims, confidence intervals and claim metadata. Both
sensitivity tables are checked cell by cell. Prediction SHA-256 fingerprints
cover every row and column, including CSV serialization. Numerical tolerances
remain rtol=1e-9, atol=1e-10. Input runtime metadata is descriptive; input bytes,
relative paths and directory membership must match exactly.

Changing these references is a change to the scientific baseline and requires
explicit review. The pipeline never rewrites them.
