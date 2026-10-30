Sparse BIC substantially reduces clean-profile over-attribution compared with MLE.
At 2,000 mutations, MLE with background selects 15.20 active signatures on average.
Sparse BIC with background selects 3.92, while SigProfilerAssignment selects 3.75.
Here, an active signature receives at least ten fitted mutations.

Removing background improves sparse BIC on these clean profiles.
Its mean attribution error reaches 0.0492, compared with 0.0496 for SigProfilerAssignment and 0.2389 for MLE without background.
The sparse-minus-SigProfilerAssignment difference is −0.00034, with a paired 95% interval of [−0.01270, 0.01331].
This comparison supports similar performance at this burden, without establishing superiority.

MLE still reconstructs observed counts more closely at 2,000 mutations.
With background, its mean cosine reaches 0.993927, compared with 0.991749 for sparse BIC and 0.990461 for SigProfilerAssignment.
Better reconstruction therefore coexists with worse recovery of the generating exposures.

At 100 mutations, SigProfilerAssignment clearly performs better on clean profiles.
Mean attribution error reaches 0.339 for SigProfilerAssignment, 0.539 for sparse BIC without background, and 0.607 with background.
At 50,000 mutations, sparse BIC without background reaches 0.00559, compared with 0.01490 for SigProfilerAssignment.
The paired difference is −0.00931, with a 95% interval of [−0.01897, −0.00126].
With background, the corresponding interval includes zero: [−0.01335, 0.00393].
These results support a burden-dependent interpretation, rather than a universal ranking.

Background helps when contamination matches its basis.
At 2,000 mutations, sparse BIC with background reaches errors of 0.0961 under uniform contamination and 0.0950 under substitution-specific contamination.
Without background, those errors reach 0.2682 and 0.3470, respectively.
For omitted signatures, background reduces sparse BIC error from 0.2866 to 0.2352.
These controls do not establish accurate recovery of arbitrary unknown processes.

Background also absorbs some known-signature mass in clean profiles.
At 2,000 mutations, sparse BIC assigns mean background mass of 0.0467.
Its clean-profile attribution error increases from 0.0492 without background to 0.0838 with background.
Background therefore needs explicit interpretation even when the solution contains few known signatures.

Across equally weighted evaluation conditions and burdens, sparse BIC with background has the lowest mean error, 0.2401.
SigProfilerAssignment reaches 0.2740, while MLE with background reaches 0.3493.
Two of four conditions deliberately match the background basis, which limits this aggregate comparison.

Several legacy objectives frequently invoke their normalized NNLS fallback.
The following rates average all evaluation conditions and burdens with background enabled.

| Method | Fits that invoke fallback |
|---|---:|
| MLE | 23.44% |
| AICc | 27.99% |
| BIC | 31.77% |
| Compat | 41.54% |
| Elastic | 95.96% |
| sparse-bic | 0% |
| sparse-aicc | 0% |

The remaining evaluated MutaGene objectives report no fallback in this panel.
All 20,736 evaluation fits return finite estimates without unresolved errors.
Fallback results remain in the comparison because they reflect current behavior.
They do not demonstrate successful optimization of the requested legacy objective.
