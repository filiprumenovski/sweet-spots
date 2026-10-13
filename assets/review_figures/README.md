# Revised review figures

These nine editorial exports are the selected 11 October figure set. They are
retained assets, with every panel mapped to its source in `manifest.json`.
`captions.json` contains the associated draft captions; S3 still needs the
original method-specific detail when integrated into the manuscript.

The production workflow regenerates the analysis figures in `results/figures/`.
It also verifies these editorial exports and their source data. A matching
source does not mean the retained artwork was freshly rendered. The reviewer
report shows both sets side by side and reports this distinction explicitly.

`recipes/` preserves historical plotting code for inspection. Machine-specific
paths are replaced with placeholders; these are not executable pipeline
entrypoints. `retained_sources/` contains the authored vector panels and Table
S3 needed to inspect the retained evidence. No provider input data are bundled.

Source hashes refer to the immutable archived results. Small CSV/JSON reference
snapshots allow a tolerance-based comparison where CPU-dependent serialization
differs while the numerical values reproduce. The two full prediction tables
are always checked byte for byte by the computational audit.
