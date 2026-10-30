A gradual transition is the strongest further-improvement candidate in these tests.
The six-method screen provides no statistically clear improvement across both low-burden metrics.
Its strongest general candidate, BIC-weighted averaging, reduces mean total error by only 0.00179 at low burden.
The descriptive interval is [-0.00526, 0.00105].
It also increases error at 100 and 999 mutations.

Fixed-background BIC helps several contaminated-profile scenarios but performs poorly on clean low-burden profiles.
For clean profiles with 100 mutations, its mean total error is 0.55917, versus the current default's 0.18319.
Signature swaps and selecting the lower-BIC initialization also produce scenario-specific regressions.
Lower model-selection scores do not guarantee lower attribution error in this screen.

The matched boundary test reveals a separate problem with the current hard threshold.
Adding one mutation increases mean exposure L1 change to 0.26723 across 64 profiles.
Keeping forward search at both counts reduces that change to 0.02014.
The largest exposure L1 change under the current rule is 1.05096.
The L1 scale ranges from zero to two for normalized weights.

The independent transition experiment improves both mean error metrics at every tested burden between 1,000 and 2,000 mutations.
All corresponding descriptive intervals exclude zero.
It preserves all 192 endpoint-control estimates exactly.
Two condition-level total-error regressions remain: 0.00083 for out-of-reference profiles at 1,001 mutations.
The other regression is 0.00109 for substitution-specific background at 1,750 mutations.
The production default remains unchanged.
The transition candidate remains experimental because it blends estimates instead of returning a likelihood-refitted selected model.
