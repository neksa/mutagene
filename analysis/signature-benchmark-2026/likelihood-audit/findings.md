### Full-model likelihood repair

All 1,536 reference fits report success without fallback.
The maximum certified remaining improvement is 0.00001564 nats per mutation.
Every archived MLE result has higher normalized negative log likelihood than the corresponding reference fit.
The smallest difference is approximately 0.000000000009 nats per mutation, which has no practical significance.

Solver success alone does not establish an accurate optimum.
Among 1,204 archived MLE fits reporting success, 900 have an optimality bound above 0.0001 nats per mutation.
The scaled-only variant exceeds that bound in 1,534 of 1,536 fits despite reporting success throughout.
A large bound alone does not prove a correspondingly large objective error.
The direct likelihood comparison and repeated-count controls supply additional evidence of suboptimal results.

### Attribution after numerical repair

The following table uses clean profiles and retains background in every MutaGene fit.

| Mutations | Estimator | Attribution error | Active known signatures | Reconstruction cosine |
|---:|---|---:|---:|---:|
| 100 | Current MLE | 0.66444 | 2.33 | 0.906760 |
| 100 | Reference full MLE | 0.66589 | 2.31 | 0.906781 |
| 2,000 | Current MLE | 0.24769 | 15.20 | 0.993927 |
| 2,000 | Reference full MLE | 0.24818 | 15.23 | 0.993921 |
| 2,000 | Sparse BIC | 0.08381 | 3.92 | 0.991749 |
| 50,000 | Current MLE | 0.05433 | 22.17 | 0.999761 |
| 50,000 | Reference full MLE | 0.05001 | 22.94 | 0.999744 |

Numerical repair improves likelihood but does not solve the main attribution problem.
At 2,000 mutations, full MLE still assigns at least ten mutations to approximately fifteen known signatures.
Sparse selection addresses this model-complexity problem separately.
Cosine can decrease when likelihood improves because these objectives differ.

### Repeated-count controls

The audit selects one clean 100-mutation profile per cancer and multiplies its counts by 1, 20, and 500.
These operations preserve channel proportions and the mathematical MLE optimum.
They do not create independent observations for scientific evaluation.

For the liver example, the original fit differs from the reference by 0.00013 in exposure L1 distance.
Multiplying every count by 500 increases that distance to 0.53357, despite a successful optimizer status.
The normalized negative log likelihood increases from 4.284886 to 4.294261.
Because the objective also worsens, nonunique equivalent mixtures do not explain this change.
The other seven profiles also show substantially worse likelihood at the largest multiplier than their reference fits.

### Parameterization ablations

These controls use one profile per condition and burden, with both background settings, yielding 24 fits per variant.
Every variant uses the same interior initialization and tight stopping tolerance.

| Objective scale | Constraint | Derivatives | Successful statuses |
|---|---|---|---:|
| Raw counts | Equality | Analytic | 8/24 |
| Per mutation | Equality | Numerical | 24/24 |
| Per mutation | Inequality | Analytic | 24/24 |
| Per mutation | Equality | Analytic | 24/24 |

These ablations use a smooth objective without the production function's conditional renormalization.
The raw-count variant often reports failure even though its normalized candidates have small optimality bounds.
Its maximum bound is approximately 0.000002 nats per mutation.
Thus, an unsuccessful status does not automatically justify discarding a candidate for NNLS.
The scaled inequality variant also converges, confirming that the inequality alone is not the central defect.
