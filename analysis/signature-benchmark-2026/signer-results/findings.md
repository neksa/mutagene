The matched benchmark completes all 216 signeR fits, but independent chains show substantial sampling instability.
These results describe this configuration and do not establish general superiority over signeR.
The comparison retains every primary fit, including fits with poor diagnostics.

### Attribution accuracy

All eight methods use the same 96 simulated profiles.
Mean total attribution error is 0.40704 for signeR and 0.33329 for signeR with background.
Joint BIC scores 0.25656, fixed BIC scores 0.21654, and SigProfilerAssignment scores 0.28403.
Full-model MLE scores 0.33982.
Lower values indicate less attribution error.

The signeR-with-background minus joint-BIC difference is 0.07673.
Its descriptive cancer-group bootstrap interval is [0.03757, 0.12310].
This interval describes the sampled profiles and does not resolve sampler instability.

The thresholded known-signature metric changes the ordering.
With that metric, signeR with background scores 0.20758, compared with joint BIC's 0.23157.
Their difference is -0.02399, with interval [-0.05547, 0.01243].
Fixed BIC scores 0.18364 and has the lowest mean under both metrics.
Two benchmark conditions directly match the fitted background family.
That design limits generalization of the fixed-background advantage.

### Mutation burden

Clean profiles with 100 mutations produce large signeR attribution errors.
Total error reaches 0.76524 without background and 0.78840 with background.
SigProfilerAssignment scores 0.24128, and joint BIC scores 0.44720.

For clean profiles with 2,000 mutations, signeR scores 0.08803 without background and 0.08724 with background.
Joint BIC scores 0.09354.
At 50,000 mutations, joint BIC scores 0.00676.
signeR scores 0.00827 without background and 0.00873 with background.
Each condition/burden cell contains only eight simulated profiles.

### Sampling limitations

The independent-chain analysis compares 24 pairs from Liver-HCC profiles.
Twenty pairs contain components above the ten-mutation activity threshold.
Nineteen of those twenty pairs have maximum active Rhat above 1.01.
Eighteen have minimum active bulk ESS below 100.
Maximum active two-chain Rhat reaches 2.00955, and minimum active bulk ESS reaches 2.73263.
Four pairs have no active components and therefore no active-component diagnostic summary.

Relative-exposure L1 differences between independent chains range from 0.00256 to 0.82578.
The median difference is 0.13966.
These discrepancies prevent a converged-posterior interpretation of the benchmark as a whole.
Longer chains and broader independent-chain checks require a separate, declared follow-up experiment.
The current comparison preserves its original settings and includes every completed fit.

### Exposure scale and runtime

Native median exposure totals differ from observed mutation burdens.
Their ratios range from 0.58414 to 0.99960 without background and 0.56353 to 1.01719 with background.
Raw-scale mean total attribution errors are 0.37748 and 0.30386, respectively.
The primary relative-exposure metrics normalize native medians by their sum.
Both scales remain available in the metrics files.

Mean fit times are 164.93 seconds without background and 171.08 seconds with background.
Different execution dates and concurrent workloads prevent a controlled runtime comparison with the preceding methods.

### Decision

This extension supplies a matched fixed-catalog signeR comparator with explicit sampling limitations.
It provides no basis for a new default change or a general superiority claim.
External-cohort validation and a calibrated context-independence null remain separate work.
The application-note manuscript still excludes this extension and the preceding sparse-selection and numerical-repair reports.
