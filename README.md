# Sweet spots

[![Release: v1.1.1](https://img.shields.io/badge/release-v1.1.1-2F6F9F)](https://github.com/filiprumenovski/sweet-spots/releases/tag/v1.1.1)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21922595.svg)](https://doi.org/10.5281/zenodo.21922595)
[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Snakemake 9.17](https://img.shields.io/badge/Snakemake-9.17-039475)](https://snakemake.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/Code-MIT-blue.svg)](LICENSE)

Code and data for *Regional sequence composition predicts O-GlcNAc organization
across species* (Rumenovski and Fehl, submitted to Cell Systems). The workflow
rebuilds the O-GlcNAc region catalog, every computational analysis and figure,
the result ledger, and an audit that checks each reported value against the
archived results.

## Reproduce everything with one command

Open a terminal on a Mac or Linux computer and paste:

```bash
git clone https://github.com/filiprumenovski/sweet-spots.git
cd sweet-spots
./reproduce
```

When the run finishes, its last line says:

```text
REPRODUCED. Every frozen scientific check passed.
```

You do not need to install Python, set up an environment or download any data:
the script fetches everything itself, checks every download, recomputes
every result from scratch, and compares each number with the published values.

| | |
|---|---|
| **You need** | macOS or Linux, an internet connection, 8 GB of RAM, 4 GB of free disk |
| **How long** | about 77 minutes on a 10-core laptop; about 20 minutes on a 96-core server (see below) |
| **If it stops** | run `./reproduce` again; it picks up where it left off |
| **If it fails** | it prints which step failed; details are in `results/logs/reproduce.log` |
| **Results** | figures in `results/figures/`, numbers in `results/analysis/`, full report in `results/review/index.html` |
| **To start over** | delete the folder and clone again |
| **Windows** | not supported directly; use WSL 2 with Ubuntu and run the same three lines there |
| **No `git`?** | download the [ZIP](https://github.com/filiprumenovski/sweet-spots/archive/refs/heads/main.zip), unzip it, open a terminal in the folder and run `bash reproduce` |

### Expected run time

| Computer | Time |
|---|---|
| Laptop, 10 cores (Apple M4) | about 77 minutes |
| Server, 96 cores (AMD EPYC) | about 20 minutes |

Fewer cores take proportionally longer. Running it again after it finishes takes a
minute or two.

<details>
<summary>What the command does, step by step</summary>

1. Checks the computer (operating system, cores, memory, disk).
2. Downloads a pinned copy of `uv` and Python 3.12.13 into `.tools/`, verifies their
   checksums, and installs every package at the exact version in `uv.lock` into
   `.venv/`. Python, conda or uv settings already on the computer are ignored.
3. Downloads all input data (the open-data release and the supplemental archive, with
   Zenodo as a backup source) and checks every file against
   `reference/input_manifest.json`.
4. Deletes the retained results and recomputes all analyses, figures and tables using
   every CPU core.
5. Runs the audit: every reported number, sensitivity row, prediction table and figure
   is compared with the archived reference. It prints `REPRODUCED` and writes
   `results/audit/PASS` only if every check passes.

Nothing is installed outside the `sweet-spots` folder.

</details>

## What is where

| Path | Contents |
|---|---|
| `workflow/`, `src/`, `config/` | Snakemake workflow, analysis code and parameters |
| `results/` | Retained analyses, tables, figures and audit (`./reproduce` recomputes them) |
| `reference/` | Archived values the audit compares against |
| `assets/review_figures/` | Designed figures; review 1–5 are paper Figures 2–6, review S1, S4, S3, S2 are paper Figures S1–S4 |
| `figure_1/` | Build for paper Figure 1 ([instructions](figure_1/README.md)) |
| `docs/input_data.md` | Input sources, versions, checksums and licenses |

For development, `make check` runs the linter, unit tests and lockfile check, and
`make dry-run` prints the full execution plan. `profiles/slurm` runs the same
workflow on a SLURM cluster.

## Licenses

The code is MIT-licensed. The supplemental input archive
([Zenodo 23330269](https://doi.org/10.5281/zenodo.23330269)) is distributed under
CC BY-NC-SA 4.0; sources and attribution are listed in
[docs/input_data.md](docs/input_data.md).

## Citation

Cite this software through `CITATION.cff` (GitHub's **Cite this repository**) or the
[concept DOI](https://doi.org/10.5281/zenodo.21922595), which resolves to the latest
release. Authors: Filip Rumenovski ([ORCID](https://orcid.org/0009-0005-1552-2826))
and Charlie Fehl ([ORCID](https://orcid.org/0000-0001-6182-0879)), Department of
Chemistry, Wayne State University.
