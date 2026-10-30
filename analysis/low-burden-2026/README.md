# Improve low-burden attribution

Adaptive joint BIC lowers low-burden error and preserves high-burden estimates in fresh synthetic validation.
The [report](report.html) contains the results, uncertainty estimates, and implementation checks.
The policy uses forward search through 1,000 mutations and retains the preceding search above that threshold.
Historical artifact hashes remain unchanged; manifest paths point to preserved source snapshots where necessary.

The sections below describe the experiment design.

## Screen candidate mechanisms

The pilot uses 64 previously examined profiles at 100 mutations.
It covers all eight cancer distributions and all four conditions.
These profiles support exploration, not independent validation.

Candidates test three mechanisms:

- AICc changes the finite-sample model penalty.
- Extended BIC adds a model-space penalty during the existing subset search.
- Sparse NNLS changes reconstruction geometry and optionally refits selected weights with multinomial MLE.

Extended BIC uses `BIC + 2 * gamma * log(binomial(p, k))`.
Here, `p` counts available components and `k` counts retained components, including background.
This exploratory mixture adaptation does not inherit a consistency guarantee from regression theory.
[Chen and Chen, 2008](https://doi.org/10.1093/biomet/asn034).

Sparse NNLS removes components when their removal increases relative Euclidean error by at most the configured tolerance.
It adds components only when the improvement exceeds twice that tolerance.
This candidate uses no tissue label or empirical exposure prior during fitting.

## Require separate validation

Select the candidate and burden policy using tuning data before generating validation results.
Evaluate both unthresholded total attribution error and thresholded known-signature error.
Report every cancer distribution and background condition, including failures.
Check intermediate burdens and high burdens separately.
Require lower low-burden error and no high-burden degradation before changing a default.
Fresh synthetic draws assess this simulation design; they do not establish patient-level generalization.
