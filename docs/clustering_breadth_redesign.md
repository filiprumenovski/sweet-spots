# Clustering-breadth redesign

## Why the original pair was retired

The former per-protein Monte Carlo count and one-dimensional “Ripley” count
were two thresholds applied to the same radius-5 close-pair statistic under the
same within-protein residue shuffle. The reported numbers were therefore not
independent evidence. The Ripley transformation was monotone in the pair count
and did not add an edge-corrected or genuinely multiscale estimand.

Per-protein significance counts were also a poor answer to the stated question.
They conflate effect breadth with power: proteins with three mapped sites have
far less chance to cross an FDR threshold than deeply sampled proteins, even
when their effect is in the same direction.

## Replacement estimands

The primary breadth scale is fixed at 10 residues, matching the manuscript's
primary self-clustering comparison. Five, 15, 20, 25 and 30 residues are
sensitivity scales and are never treated as independent tests.

The residue-opportunity analysis reports:

1. the protein-population fraction with observed close-pair fraction greater
   than its exact fixed-count expectation;
2. among proteins carrying both marks, the fraction in which O-GlcNAc excess is
   greater than phosphorylation excess; and
3. the fraction satisfying each directional criterion at every sensitivity
   scale.

Wilson intervals describe protein-population fractions. Whole-protein
bootstrap intervals remain the uncertainty measure for ratios of mean observed
to mean expected pair fractions.

## Peptide-exposure-conditioned null

The detection-aware sensitivity analysis was specified as follows before its
new implementation was run:

- use the strict sequence-validated O-GlcNAc atlas;
- require at least three unique sites per protein;
- use 10 residues as the primary radius;
- map every reported peptide exactly to the canonical protein sequence;
- use the deterministic theoretical tryptic interval only when no reported
  peptide maps to a site;
- partition each protein's serines and threonines by residue identity and the
  exact set of observed peptide intervals covering that residue;
- redraw the observed number of sites independently within each disjoint
  exposure stratum, without replacement;
- preserve every protein's total site count exactly in every simulation;
- use 2,000 simulations per protein and deterministic accession-derived seeds;
- report the ratio of equal-protein mean observed and null pair fractions,
  alongside the protein-population fraction with positive excess;
- treat all non-primary radii and any region-count endpoint as sensitivity
  analyses.

This null is intentionally conservative. A site in a stratum with no unused
eligible residue is fixed rather than moved to a chemically or detectably
different location. The analysis therefore asks whether clustering remains
after conditioning on the peptide exposure actually present in the catalogue,
not whether every site could have been observed everywhere in the protein.

## Independent-evidence restrictions

Two complementary robustness populations address the more direct catalogue
artifact concern: that clustering could be carried by results from one study or
one repeatedly reported peptide. They retain only sites supported by at least
two distinct PubMed records or by at least two distinct, exactly mapped peptide
sequences, respectively. Both use the same exact fixed-count all-serine/threonine
expectation, the same fixed primary radius and the same whole-protein bootstrap
as the primary analysis.

These restrictions are expected to look cleaner if the phenomenon is genuine:
they trade catalogue breadth for better-supported site localization. They are
post hoc robustness analyses, however. Their effect sizes cannot be compared
with the unrestricted cohort as though the populations were interchangeable.

## Guardrails

- The redesign is labeled post hoc and cannot be represented as prespecified.
- No null, radius or population is selected according to the resulting effect
  size.
- The old Monte Carlo and pseudo-Ripley counts must not be described as
  independent confirmation.
- The detection-aware result is a sensitivity analysis for ascertainment bias,
  not proof of a biochemical clustering mechanism.
- Retained output must pass exact site-count, unique-position, shard-coverage
  and receipt-digest checks before entering the result ledger.
