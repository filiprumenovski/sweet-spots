# Region definition

The paper uses *consensus region* as the name of a fixed, operational two-stage
object. A region is a variable-width interval spanned by a locally dense chain
of validated O-GlcNAc sites. It is not a fixed-width window, a predicted domain,
or a claim about a physical boundary.

## Stage 1: strict core membership

For each protein, the workflow sorts the unique, sequence-validated O-GlcNAc
serine and threonine coordinates. Single-linkage segmentation joins consecutive
sites when their coordinate difference is no greater than the maximum gap. A
component qualifies only when it contains at least three sites.

The archived implementation runs this operation at gaps of 5, 8, 10, 12, and
15 residues, then retains sites that belong to a qualifying component at every
gap. These thresholds are nested: increasing the gap can only merge or extend a
component. The intersection is therefore exactly the set selected at the
strictest five-residue gap. We call these retained coordinates *strict core
sites*.

The five configured passes are preserved for exact provenance, but they are not
independent votes and do not make core membership parameter-free.

## Stage 2: reported region grouping

The workflow re-segments the strict core sites with a ten-residue maximum gap,
again requiring at least three sites. Each resulting component becomes one
reported region. Its `start` and `end` are the first and last core-site
coordinates, both one-based and inclusive. Its `span` is `end - start + 1`, and
its `valence` is the number of core sites in the component.

This second pass can merge strict core chains separated by 6 to 10 residues.
It adds no padding around the endpoints. Regions have variable widths, and
single-linkage grouping imposes no maximum length.

## Worked example

Consider validated sites at positions 100, 103, 108, 116, 120, and 124.

- At gap 5, the sites form two qualifying strict core chains: 100-103-108 and
  116-120-124.
- At the final gap 10, the eight-residue break between positions 108 and 116 is
  bridged.
- The reported region is therefore the inclusive interval 100-124, with span 25
  and valence 6.

## Realised catalogue and sensitivity

The retained catalogue contains 824 regions on 477 proteins. It contains 4,196
strict core sites from 12,491 validated human O-GlcNAc sites. Figure S4 reports
a post hoc factorial sensitivity grid with strict-core gaps of 3 through 7,
minimum component sizes of 3 through 5, and final grouping gaps of 8, 10, and
12. The complete 45-definition catalogue grid is retained in
`results/analysis/region_definition_sensitivity/catalogue_grid.csv`.

The composition-only, covariate-only, and adjusted-composition models are
refitted for the 15 strict-core definitions at final gap 10. Cross-validation
holds out entire proteins, and uncertainty in the adjusted composition
increment comes from paired whole-protein bootstrap resampling. Those model
results are retained in
`results/analysis/region_definition_sensitivity/model_grid.csv`. Because this
grid was designed after inspecting the primary analysis, it is labelled as a
post hoc robustness analysis rather than an independent confirmation.

The executable definition is
[`consensus_regions.py`](../src/regional_code_paper/analysis/consensus_regions.py),
and the frozen identifier is
`og_consensus_gaps5-8-10-12-15_core-gap10_min3_strict-v1`.
