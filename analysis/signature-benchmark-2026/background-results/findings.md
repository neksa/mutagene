Sparse BIC with joint background selection now supplies the shared Python, CLI, and web default.
It uses constrained MLE for every candidate's weights.
It counts retained known and background components in BIC.
Users can retain full-model fitting with `--method MLE` or `func="MLE"`.
Users can retain the previous sparse background policy with `--background-selection fixed`.

### Why this default

Joint selection lowers total attribution error against full-model MLE in all twelve condition/burden cells of both panels.
The fresh balanced-panel mean decreases from **0.35626** to **0.26616**.
This comparison uses unthresholded known-signature and aggregate background-mass error.
It retains the preceding reassessment's total-attribution metric.

For fresh clean profiles with 2,000 mutations, total attribution error decreases from **0.25692** to **0.05823**.
The joint-minus-MLE difference is **−0.19869**, with a within-cancer bootstrap interval of **[−0.22052, −0.17623]**.
Joint selection retains 4.50 active known signatures on average, compared with MLE's 15.22.
Reconstruction cosine decreases from 0.99327 to 0.99035 as the model becomes sparser.

Background selection addresses the specific concern about excess background attribution.
At 2,000 mutations, fresh clean background mass decreases from **0.05444** under fixed BIC to **0.00168** under joint BIC.
The paired panel shows the same direction: **0.04565** decreases to **0.00022**.
The true background mass is zero in these clean simulations.
Reference-version mismatch and sampling noise still affect fitted attribution.

The default favors explicit background selection and clean-profile interpretation over a fixed six-component nuisance model.
It does **not** minimize every benchmark metric or dominate every alternative.

### Tradeoffs and limits

Fixed-background BIC achieves the lowest balanced-panel total attribution error: **0.24833**, compared with joint BIC's **0.26616**.
Two of four conditions directly match the fixed background family.
That design favors background modeling, but the fixed-policy advantage remains a real measured result.
For fresh uniform-background profiles at 2,000 mutations, fixed BIC scores **0.11305** and joint BIC scores **0.17026**.
For substitution-specific background at 50,000 mutations, fixed BIC scores **0.01289** and joint BIC scores **0.02663**.
Joint selection can remove weak background components and underfit genuine diffuse background.
The fixed policy remains available when that tradeoff matters.

SigProfilerAssignment remains substantially stronger at low burden.
For fresh clean profiles with 100 mutations, its total attribution error is **0.28875**, versus joint BIC's **0.50179**.
Joint BIC and SigProfilerAssignment have similar clean-profile accuracy at 2,000 mutations.
Their total-error difference is **0.00195**, with interval **[−0.01930, 0.02153]**.
Joint BIC performs better at 50,000 mutations in the fresh panel.
That advantage does not clearly reproduce in the paired panel.
These results do not establish universal superiority over SigProfilerAssignment.

Thresholded known-signature error gives a different low-burden comparison.
For fresh clean profiles with 100 mutations, MLE scores **0.42486**, while joint BIC scores **0.45346**.
Joint BIC also scores worse than MLE in the other three 100-mutation cells under this metric.
The ten-mutation threshold removes many small MLE exposures from this comparison.
Total attribution error includes those unthresholded exposures and background-mass error.
The report supplies both metrics instead of treating them as interchangeable.

Block selection removes all background in every clean profile in these panels.
Its fresh balanced-panel total attribution error is **0.28617**, worse than joint selection's **0.26616**.
Joint selection offers finer background control and lower computation cost than the block comparison here.

Mean fresh-panel fitting times are 2.06 seconds for joint BIC and 0.08 seconds for MLE.
Block selection averages 3.48 seconds.
Concurrent execution limits hardware-independent runtime comparisons.
Bootstrap analysis repeats subset selection for each replicate and therefore costs more than full-model MLE.

### Numerical verification

Both panels finish **4,608 fits**, totaling **9,216 configuration fits**, with no reported errors.
All MutaGene configurations report zero fallbacks.
Independent checks validate **7,680 final-model likelihood certificates** and **6,144 sparse-model BIC calculations**.
The maximum recorded refit gap is **1.9983e-5**, below the **2e-5** tolerance.
The verifier also recalculates 400 metric records and confirms exact paired input identity.
Repeated MLE, fixed BIC, no-background BIC, and SigProfilerAssignment attribution scores match the preceding panel.

The solver certificate bounds likelihood optimization error within each selected model.
It does not certify global subset selection, biological correctness, or exposure uncertainty.
The default decision uses simulated profiles and requires separate external-cohort validation for generalization claims.

### Regression and package checks

The full suite passes **513 tests** with **68.97%** coverage.
This change adds sixteen cases beyond the preceding likelihood reassessment.
Tests cover background removal, background-only fits, BIC dimensions, block comparisons, unsupported models, determinism, and mutation-burden effects.
Interface tests check CLI/Python agreement and pooled and per-sample web defaults.
Tests retain explicit failure propagation and full-model MLE count-scaling checks.
Twelve current-default fits match benchmark estimates across all condition/burden cells.

Strict Sphinx documentation checks and Quarto report rendering pass.
A clean installed wheel passes the demo, analytic sparse fit, web home-page, and static-asset checks.
Ruff, Black, and whitespace checks pass.
The report retains previous evidence through source snapshots and hash-preserving manifest relocation.
