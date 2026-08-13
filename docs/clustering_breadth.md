# Clustering breadth

This analysis asks whether regional O-GlcNAc self-clustering is broadly
distributed across proteins and robust to peptide-level detection opportunity.

## Estimands

The primary scale is fixed at 10 residues, matching the primary self-clustering
comparison. Radii of 5, 15, 20, 25, and 30 residues are sensitivity analyses,
not independent tests.

The residue-opportunity analysis reports:

1. the fraction of proteins whose observed close-pair fraction exceeds its
   exact fixed-count expectation;
2. among proteins carrying both marks, the fraction in which O-GlcNAc excess
   exceeds phosphorylation excess; and
3. the fraction satisfying each directional criterion at every sensitivity
   scale.

Wilson intervals describe protein-population fractions. Whole-protein
bootstrap intervals quantify uncertainty for ratios of mean observed to mean
expected pair fractions.

These endpoints replace per-protein significance counts, which conflate effect
breadth with statistical power. A protein with three mapped sites has much less
power to cross an FDR threshold than a deeply sampled protein even when the
effect points in the same direction.

## Peptide-exposure-conditioned null

The detection-aware sensitivity analysis follows a fixed procedure:

- use the strict, sequence-validated O-GlcNAc atlas;
- require at least three unique sites per protein;
- use 10 residues as the primary radius;
- map every reported peptide exactly to the canonical protein sequence;
- use the deterministic theoretical tryptic interval only when no reported
  peptide maps to a site;
- partition each protein's serines and threonines by residue identity and the
  exact set of observed peptide intervals covering that residue;
- redraw the observed number of sites independently within each disjoint
  exposure stratum, without replacement;
- preserve every protein's site count in every simulation;
- use 2,000 simulations per protein and deterministic accession-derived seeds;
- report both the equal-protein observed-to-null ratio and the fraction of
  proteins with positive excess; and
- treat all non-primary radii and region-count endpoints as sensitivities.

The null is deliberately conservative. A site in a stratum with no unused
eligible residue remains fixed rather than moving to a chemically or detectably
different location. The test therefore asks whether clustering remains after
conditioning on the peptide exposure represented in the catalogue.

## Independent-evidence restrictions

Two robustness populations address repeated reporting from a single study or
peptide. They retain sites supported by at least two distinct PubMed records or
by at least two distinct, exactly mapped peptide sequences. Both use the same
fixed-count serine/threonine expectation, primary radius, and whole-protein
bootstrap as the primary analysis.

These restrictions exchange catalogue breadth for stronger localization
support. They are post hoc robustness analyses, so their effect sizes are not
compared with the unrestricted cohort as if the populations were identical.

## Interpretation guardrails

- The robustness extensions are labeled post hoc, not prespecified.
- No null, radius, or population is selected according to the resulting effect
  size.
- Closely related pair-count transformations are not presented as independent
  confirmation.
- The detection-aware result addresses ascertainment bias; it does not by
  itself establish a biochemical mechanism.
- Outputs must pass exact site-count, unique-position, shard-coverage, and
  receipt-digest checks before entering the result ledger.
