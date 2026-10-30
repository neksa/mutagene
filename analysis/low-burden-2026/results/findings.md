Adaptive search improves low-burden attribution while preserving every high-burden estimate in the validation panel.
At 100 mutations, mean total attribution error decreases from 0.55723 to 0.38214, a 31.4% reduction.
The reductions are 28.0% at 250 mutations, 24.7% at 500 mutations, and 20.8% at 1,000 mutations.
Thresholded known-signature error also decreases at every tested low burden.

The frozen policy passes all eight low-burden confidence-interval checks.
At 100 mutations, the paired total-error difference is -0.17510, with interval [-0.21647, -0.13481].
Both mean errors improve in all sixteen low-burden condition/burden cells.
These results support the rule within this synthetic design.

All 512 profiles above 1,000 mutations retain exactly identical fitted weights.
The default uses the preceding search above that threshold.
This preserves high-burden behavior by construction, rather than relying on an average error comparison alone.

The sensitivity panel removes the generating-truth cutoff.
Total error decreases from 0.67708 to 0.51561 at 50 mutations.
At 100 mutations, it decreases from 0.62459 to 0.42439.
Both sensitivity burdens also pass the paired interval checks for both metrics.

The same-run timing comparison favors forward search.
At 100 mutations, mean fitting time decreases from 0.486 seconds to 0.172 seconds.
At 1,000 mutations, it decreases from 1.116 seconds to 0.306 seconds.
Concurrent workloads limit interpretation of these runtime measurements.

This result improves the preceding MutaGene default on fresh synthetic draws.
It does not establish patient-level accuracy or superiority over external tools.
