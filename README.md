# Sweet spots: reproducibility code

This repository is the paper-only computational record for *Sweet spots: how a
promiscuous glycosyltransferase reads a conserved regional code*. It rebuilds
the consensus region object, every load-bearing analysis retained in the Cell
Systems manuscript, all main and supplemental figures, a tidy result ledger,
and a final numerical and figure-integrity audit.

Exploratory, retired, disease, network, LLPS, stoichiometry, and legacy cassette
analyses from the parent research repository are intentionally absent.

## Input data

The workflow starts from a frozen, analysis-ready data bundle. It does not
download data or call a remote service. `source_root` in `config/config.yaml`
must point to a directory with the following files:

| Path below `source_root` | Contents | Used for |
|---|---|---|
| `data/interim/fasta_human.parquet` | Canonical human protein sequence spine | Consensus regions, self-clustering, regional models |
| `data/interim/fasta_all.parquet` | Canonical multispecies protein sequence spine | Scanner validation, evolutionary transfer, FG-NUP background |
| `data/interim/iupred_residue_scores.parquet` | Residue-level IUPred disorder scores | Covariate control, scanner features, matched FG-NUP analysis |
| `analysis/revalidation/data/atlas_unambiguous.csv` | Strict sequence-validated O-GlcNAc site evidence with species and publication support | Region definition, labels, evolutionary comparisons |
| `data/processed/landscape/ptm_unified.parquet` | Canonical O-GlcNAc, phosphorylation, acetylation, and ubiquitination site catalogue | Self-clustering and O-GlcNAc/phosphorylation analyses |
| `data/processed/regions/oglcnac_consensus_regions.parquet` | Archived canonical consensus-region object | Independent checksum comparison of the rebuilt object |
| `data/interim/msa/*.afa` | Per-protein human-to-ortholog alignments | Alignment-based site and regional conservation |
| `data/external/multispecies_oglcnac/rice_sequences.fasta` | Rice sequences linked to observed O-GlcNAc sites | Alignment-free cross-species transfer |
| `data/external/multispecies_oglcnac/arabidopsis_sequences.fasta` | Arabidopsis sequences linked to observed O-GlcNAc sites | Alignment-free cross-species transfer |
| `data/external/multispecies_oglcnac/drosophila_sequences.fasta` | Drosophila sequences linked to observed O-GlcNAc sites | Alignment-free cross-species transfer |
| `data/external/multispecies_oglcnac/celegans_sequences.fasta` | *C. elegans* sequences linked to observed O-GlcNAc sites | Alignment-free cross-species transfer |
| `analysis/revalidation/data/ogt_orthologs/O15294.fasta` | Human OGT sequence | Mapping structural residues to the canonical sequence |
| `analysis/revalidation/data/ogt_orthologs/orthologs_aln.fasta` | OGT ortholog alignment | Conservation of substrate-reading and catalytic surfaces |
| `analysis/revalidation/data/ogt_structures/*.cif` | Experimentally determined OGT structures | Peptide-contact, surface, core, and interface assignments |

The frozen development bundle contains 16,450 per-protein ortholog alignments
and 44 OGT structure files. Kinase specificity matrices are supplied by the
pinned `kinase-library==1.7.1` dependency rather than by a separate data file.
The first workflow rule verifies that every required path exists and writes
relative paths, byte counts, SHA-256 hashes, and directory-tree hashes to
`results/provenance/input_manifest.json`. This manifest is the authoritative
identity record for the frozen bundle. The ambiguous-site atlas, raw downloads,
and all retired intermediate products are deliberately outside the input
contract.

Validate and fingerprint the input bundle without launching an analysis:

```bash
uv run snakemake --snakefile workflow/Snakefile --profile profiles/local \
  results/provenance/input_manifest.json
```

Expected layout:

```text
source_root/
├── analysis/revalidation/data/
│   ├── atlas_unambiguous.csv
│   ├── ogt_orthologs/{O15294.fasta,orthologs_aln.fasta}
│   └── ogt_structures/*.cif
└── data/
    ├── external/multispecies_oglcnac/*_sequences.fasta
    ├── interim/{fasta_human.parquet,fasta_all.parquet,iupred_residue_scores.parquet}
    ├── interim/msa/*.afa
    ├── processed/landscape/ptm_unified.parquet
    └── processed/regions/oglcnac_consensus_regions.parquet
```

## Run

The workflow expects the frozen analysis-ready data tree named by
`config/config.yaml`. By default, it looks under `inputs/`. For an archive or a
reviewer machine, either unpack the bundle there or change only `source_root`
to the unpacked data directory.

```bash
cd sweet-spots
uv sync --frozen
uv run snakemake --snakefile workflow/Snakefile --profile profiles/local
```

The local profile defaults to four concurrent cores and a 32 GB aggregate
memory budget. Override either explicitly for a larger workstation, for example
`--cores 16 --resources mem_mb=96000`. The cluster profile can submit up to 128
independent jobs and lets SLURM enforce each rule's declared resources.

On a SLURM cluster, set any site-specific account/partition flags on the command
line or in a private profile layered over the included portable defaults:

```bash
uv run snakemake --snakefile workflow/Snakefile --profile profiles/slurm
```

For a clean-room rerun:

```bash
uv run snakemake --snakefile workflow/Snakefile --profile profiles/local --delete-all-output
uv run snakemake --snakefile workflow/Snakefile --profile profiles/local
```

Snakemake treats the source tree as immutable. All generated files live under
`results/`. The principal products are:

- `results/tables/manuscript_result_ledger.parquet`: one row per manuscript-facing
  estimand, materialised with the DuckDB CLI and Zstandard compression.
- `results/figures/`: vector PDF and 400 dpi PNG for five main and three
  supplemental figures.
- `results/audit/reproduction_report.json`: numerical and figure-integrity
  checks.
- `results/provenance/input_manifest.json`: SHA-256 hashes of every
  analysis-ready input.

## Design boundaries

Every stochastic split is grouped by protein. Model selection and score
threshold calibration for the regional scanner occur within outer-training
partitions. The held-out outer fold is used once. Bootstrap and permutation
resampling preserve proteins as the inferential unit unless a result is
explicitly structural or catalogue-level.

Python modules exchange transparent CSV and JSON products. Parquet
materialisation is performed only by the DuckDB CLI rule, keeping the binary
storage operation declarative and reviewable. Figure code consumes only files
declared in the Snakemake DAG and never types a result from the manuscript.

## Repository layout

```text
sweet-spots/
├── config/                   # validated scientific and execution parameters
├── docs/                     # architecture and failure guarantees
├── profiles/                 # portable local and SLURM execution defaults
├── results/                  # retained reviewer-facing artifacts only
├── src/regional_code_paper/
│   ├── analysis/             # scientific estimands
│   ├── core/                 # configuration, randomness, atomic I/O
│   ├── execution/            # restartable map/reduce commands
│   ├── models/               # feature extraction and predictive models
│   ├── reporting/            # manifests, figures, ledgers, audits
│   └── sql/                  # declarative DuckDB analyses
├── tests/unit/               # tests mirror the source architecture
└── workflow/                 # Snakemake DAG and domain-specific rules
```

## Parallel architecture

The domain-specific rule files in `workflow/rules/` mirror the paper's dependency
structure. Expensive stages use explicit map/reduce boundaries:

- self-clustering maps independently over the four PTMs;
- the peptide-exposure-conditioned clustering null maps proteins over
  configurable shards;
- regional prediction maps over five protein-grouped outer folds;
- kinase-library scoring maps matched peptide pairs over configurable shards;
- scanner validation maps five outer folds, while keeping feature extraction in
  one shared, immutable cache.

Shard counts live in `config/config.yaml`. They affect scheduling only. Record
ownership uses a sorted global index, and random seeds use that same global
index, so a 4-shard and a 64-shard run target the same numerical result.

Every map job commits outputs atomically and then writes a receipt containing
row count, byte count, and SHA-256. Reducers fail closed on missing receipts,
duplicate keys, schema drift, or incomplete fold coverage. Each expensive job
also has a separate log and Snakemake benchmark file. Interrupted jobs can be
resubmitted safely with `--rerun-incomplete`; completed shards are not repeated.
The complete map/reduce contract and failure model are documented in
[`docs/architecture.md`](docs/architecture.md).
The rationale, fixed endpoints and reporting guardrails for the clustering
breadth redesign are documented in
[`docs/clustering_breadth_redesign.md`](docs/clustering_breadth_redesign.md).

Before spending compute, inspect the resolved graph:

```bash
uv run snakemake --snakefile workflow/Snakefile --profile profiles/local --dry-run
uv run snakemake --snakefile workflow/Snakefile --rulegraph | dot -Tpdf > dag.pdf
```

Quality gates are intentionally ordinary:

```bash
uv run --group dev ruff check src tests
uv run --group dev pytest -q
uv lock --check
uv run snakemake --snakefile workflow/Snakefile --dry-run --forceall --cores 32
```

Snakemake's generic linter recommends a per-rule Conda environment. This repo
deliberately uses one exact, checked-in `uv.lock` instead: every rule invokes
`uv run --frozen`, so dependency drift fails rather than being resolved during
a job. The SLURM executor plugin is pinned in that same lock.

## Audit-sensitive distinctions

The consensus object contains 824 regions on 477 region-bearing proteins;
3,535 is the validated scanner universe. The strict self-clustering cohort
contains 12,491 validated sites, whereas 12,501 belongs to a different scanner
input universe. The final audit fails these common unit substitutions.

The archived self-clustering resampling uses seed 20260814. Other stochastic
analyses derive deterministic analysis-specific seeds from base seed 42. Both
values are exposed in `config/config.yaml` and in result provenance.
