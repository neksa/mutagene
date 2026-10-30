# Explore extended-context logos

Status: feature proposal; the current implementation does not provide this display.

## Scientific question

Do mutations associated with an SBS96 signature show additional sequence preferences beyond its immediate flanking bases?
Display five-base contexts for ±2 flanks and seven-base contexts for ±3 flanks.
Keep the central substitution explicit rather than hiding it inside a reference-sequence logo.

## Proposed calculation

Use each sample's fitted signature exposures and SBS96 channel probabilities to calculate fractional mutation responsibilities.
Include the context-independent background components in the responsibility denominator.
Retain background responsibility instead of redistributing it across named signatures.
For each signature, aggregate extended reference contexts using those fractional responsibilities as weights.
Orient each mutation consistently around a pyrimidine reference base.
Use the same assembly, callable regions, and sequence filters as the profile analysis.

A labelled count matrix alone cannot reconstruct these flanks.
The feature requires mutation coordinates and a reference genome.
Record coordinate conventions, strand transformations, assembly, and excluded contexts.

## Interpretation and controls

The central trinucleotide already drives signature assignment.
Its apparent enrichment therefore does not independently validate the assignment.
Outer flanks add information, but genomic composition and correlated signature assignments can still confound the result.

Show both observed base frequencies and enrichment against matched genomic opportunities.
Match background opportunities on the central trinucleotide and relevant callable regions.
Keep substitution classes separate when constructing matched backgrounds.
Compare coding and neutral regions separately when selection can alter motif frequencies.
Test sensitivity to recurrent driver hotspots.
Report missing sequence, raw mutation counts, total fractional weight, and assignment uncertainty.
Use sample-level resampling to preserve within-sample dependence.
Flag profiles with insufficient effective support.

The display describes mutations associated with a signature under the fitted model.
It does not establish a new SBS1536 or SBS24576 signature.
It does not prove a mutational mechanism.

## Initial acceptance checks

- Recover a simulated outer-flank motif with sufficient support.
- Show no outer-flank enrichment under a matched composition-only control.
- Preserve reverse-complement equivalence and assembly provenance.
- Retain background assignments and report ambiguous signature responsibility.
- Distinguish five-base and seven-base contexts without changing the original SBS96 fit.

## Delivery boundary

Implement this feature separately from sparse signature selection.
Use Vega-Lite for interactive summaries where appropriate.
Preserve tabular counts and opportunity denominators for independent reproduction.
