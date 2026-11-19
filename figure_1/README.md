# Figure 1 (concept figure)

Figure 1 is assembled outside the main workflow, so `./reproduce` does not build it. Everything the
final assembly reads is in this folder. Running

```bash
python3 build_figure_1_concept.py
```

writes `Figure_1_Concept_v2.svg`, `.pdf` and `.png`. It needs Pillow and `rsvg-convert`. The SVG
is byte-identical to the submitted Figure 1.

## Panels and their inputs

| Panel | Built by | Inputs |
|---|---|---|
| A, OGT dimer | `render_ogt_assembly.py` (PyMOL 3.1.0) | `structures/7YEA.cif`; `structures/3PE4.cif` superposed for UDP and acceptor peptide |
| A, reaction | `render_chemistry.py` (RDKit 2026.09.1) | none |
| B, substrate map | `compute_substrate_map.py` → `data/substrate_map.json` | O-GlcNAc Database (OVSlab) all-species export, accessed 4 June 2026 (`OVSLAB_CSV`); repository inputs; `data/uniprot_human_subcellular.tsv` (query in `data/uniprot_human_subcellular.query.txt`) |
| C, FOXK1 model | `render_foxk1_model.py`, `render_glcnac_site.py` (PyMOL) | `structures/AF-P85037-F1-model_v6.pdb`; `structures/5UN8.pdb` (GlcNAc-Ser unit) |
| C, residue track | `compute_track.py P85037 FOXK1` → `data/foxk1_track.json` | PhosphoSitePlus `Phosphorylation_site_dataset.gz` (`PSP_SITES`); repository inputs |
| D, transfer | values in `build_figure_1_concept.py` | result ledger rows `transfer.human_to_target.rice` and `transfer.rice_to_human.human` (Figure 4A, Table S5) |

"Repository inputs" are the frozen files that `./reproduce` downloads into `inputs/`:
- `data/interim/fasta_human.parquet`
- `analysis/revalidation/data/atlas_unambiguous.csv`
- `data/processed/regions/oglcnac_consensus_regions.parquet`

The O-GlcNAc Database export and PhosphoSitePlus are provider-restricted, so they are not
redistributed here. Set `OVSLAB_CSV` and `PSP_SITES` to your own copies to rerun the two `compute_*`
steps. Their outputs are committed under `data/`, so the final assembly runs without them.

The human and rice drawings in `assets/organisms/` are carried over from the earlier Figure 1
assembly. `Assembly_Record.json` records the canvas, the sources and the transfer estimates of the
last build.
