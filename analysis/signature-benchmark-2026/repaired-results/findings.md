### Numerical reliability improves

Both panels finish without unresolved errors: 20,736 paired fits and 22,272 fresh fits.
Every repaired MLE and sparse fit completes without an NNLS fallback.
Independent verification checks all 9,216 full-MLE configuration certificates across both panels, including the compatibility aliases.
It also checks all 6,144 sparse fits against their recorded maximum refit bounds.
The largest full-MLE bound is 0.00001564 nats per mutation.
The largest sparse-refit bound is 0.00001999, below the declared 0.00002 tolerance.

At 50,000 mutations with background, pre-repair MLE invokes fallback in 69.53% of paired fits and 68.75% of fresh fits.
Repaired MLE invokes fallback in zero fits in both panels.
Six paired sparse-BIC fits and thirteen fresh sparse-BIC fits require one or more subset-refit retries.
All finish successfully within the likelihood model.

The original count-scaling regression now produces identical weights at multipliers of 1, 20, and 500.
Its normalized negative log likelihood remains 4.284886034553 throughout.
The paired SigProfilerAssignment rerun also reproduces the earlier attribution errors exactly.

### Attribution gains from numerical repair remain modest

For fresh clean profiles with 2,000 mutations and background, MLE error decreases from 0.24049 to 0.23824.
The repaired-minus-original difference is −0.00225, with a paired 95% interval of [−0.00557, −0.00018].
At 50,000 mutations, error decreases from 0.05649 to 0.05031.
That difference is −0.00618, with an interval of [−0.00898, −0.00350].
At 100 mutations, the interval includes zero.

Numerical repair therefore corrects the estimator without resolving the main over-attribution problem.
Fresh clean profiles with 2,000 mutations still receive approximately fifteen active signatures under full MLE.
Sparse BIC with background selects 3.94 active signatures and reaches error 0.08339.
Without background, sparse BIC selects 4.25 and reaches error 0.04649.

### The fresh panel limits the earlier superiority claim

SigProfilerAssignment reaches error 0.05121 on fresh clean profiles with 2,000 mutations.
Sparse BIC without background differs by −0.00472, with a paired 95% interval of [−0.01855, 0.00856].
This interval supports similar performance without establishing superiority.

At 50,000 mutations, background-free sparse BIC reaches 0.00661, while SigProfilerAssignment reaches 0.00756.
Their difference is −0.00095, with an interval of [−0.00228, 0.00024].
The earlier paired high-burden advantage therefore does not clearly generalize to these fresh draws.
With background, sparse BIC reaches error 0.01365 and performs worse than SigProfilerAssignment on these clean high-burden profiles.

At 100 mutations, SigProfilerAssignment remains substantially better.
Its mean error is 0.33306, compared with 0.49488 for sparse BIC without background and 0.55081 with background.

### Background has a condition-dependent benefit

For fresh profiles with 2,000 mutations, background improves sparse BIC under both matched contamination controls.

| Condition | Sparse BIC without background | Sparse BIC with background | SigProfilerAssignment |
|---|---:|---:|---:|
| Clean | 0.04649 | 0.08339 | 0.05121 |
| 20% omitted-signature contribution | 0.26696 | 0.22649 | 0.28164 |
| 20% uniform contamination | 0.26644 | 0.10142 | 0.23032 |
| 20% substitution-specific contamination | 0.33430 | 0.08254 | 0.32458 |

In clean profiles, sparse BIC assigns mean background mass of 0.04565 despite zero background truth.
This explains part of its larger attribution error when background remains available.
Across the balanced fresh panel, sparse BIC with background averages error 0.24412, compared with 0.28382 for SigProfilerAssignment.
Two of four conditions deliberately match the background basis, limiting any aggregate ranking.

### Legacy objectives still require caution

The fresh panel retains normalized-NNLS fallback in several non-MLE objectives with background enabled.

| Configuration | Fallback rate |
|---|---:|
| Pre-repair MLE | 23.31% |
| Repaired MLE, MLE-scaled, and Compat | 0% |
| Sparse BIC and sparse AICc | 0% |
| Legacy AICc | 25.52% |
| Legacy BIC | 34.11% |
| Elastic | 17.06% |

The remaining evaluated objectives report no fallback in this fresh panel.
These legacy fallback results remain visible because they describe actual behavior.
They do not establish successful optimization of the requested legacy objective.

The evidence supports repaired MLE as a numerically reliable implementation of its stated full-model likelihood.
It supports sparse selection as a separate attribution improvement under the tested conditions.
It does not support a universal default switch to sparse BIC across mutation burdens and contamination settings.
