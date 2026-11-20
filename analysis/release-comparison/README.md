# Compare the release default with other attribution methods

## Paired signatures from the fitting catalog

The pilot compares 32 simulated SBS96 profiles at two mutation burdens.
Seed 20261120 selects 16 signature pairs from the bundled COSMICv3 reference.
Each mixture assigns 70% to the first signature and 30% to the second.
Each pair generates one multinomial draw at each burden.
Both methods receive identical counts and all 67 reference columns.
The generator contains no additional background or missing signatures.

| Mutations per sample | Profiles | MutaGene mean error | SigProfilerAssignment mean error |
|---:|---:|---:|---:|
| 100 | 16 | 0.08452 | 0.22125 |
| 2,000 | 16 | 0.01677 | 0.01863 |

Error equals half the summed absolute differences between estimated and generating weights, including aggregate background mass.
Lower values indicate closer agreement.
MutaGene uses its default joint-background sparse BIC fit.
SigProfilerAssignment 1.1.5 uses its default assignment settings with the supplied reference.
The adapter aligns both input tables by channel label.

MutaGene has lower error in 12 of 16 low-burden pairs. One pair ties.
At 2,000 mutations, nine pairs favor MutaGene and seven favor SigProfilerAssignment.
The generator and both fitters share a catalog.
The experiment tests clean two-component mixtures only.

## Empirical exposure mixtures across mutation burdens

The archived panel contains 3,072 simulated profiles from eight cancer exposure distributions.
Each of 48 condition and burden combinations contains 64 profiles.
The generator uses COSMIC v3.3.1 and preserves all generating exposures.
Each method fits the same counts against COSMIC v3.
Two of the four simulation conditions match MutaGene's background components.

| Condition | Mutations | MutaGene default | SigProfilerAssignment | Difference (95% paired interval) |
|---|---:|---:|---:|---:|
| All conditions and burdens | 10 to 10,000 | 0.32691 | 0.36866 | −0.04174 |
| Clean | 10 | 0.76321 | 0.66639 | +0.09682 (+0.03834 to +0.15461) |
| Clean | 100 | 0.3116 | 0.4033 | −0.09164 (−0.14843 to −0.03574) |
| Clean | 2,000 | 0.0605 | 0.0584 | +0.00209 (−0.01685 to +0.01979) |
| 20% substitution-specific background | 2,000 | 0.1467 | 0.2908 | −0.14412 (−0.17227 to −0.11662) |

Differences subtract SigProfilerAssignment error from MutaGene error.
Negative differences favor MutaGene.
Of 48 paired intervals, 23 favor MutaGene, three favor SigProfilerAssignment, and 22 include zero.
All three intervals favoring SigProfilerAssignment occur at ten mutations.
The intervals use 2,000 within-cancer bootstrap resamples without multiple-comparison adjustment.

The [complete table](burden-panel-means.tsv) includes full-model MLE and full-start BIC as internal baselines.
We recalculated all means from 12,288 per-fit records.
The archived hashes for `sparse.py` and `likelihood.py` match the release candidate.
Local Git object `b4dd06407bc61d985c9acee69cfe04d39b0dde5f` contains the original `analysis/burden-range-2026/results` directory.
Its `metrics.tsv` has SHA256 `25e1f2fa678fa94e8b1850d322335767cc2198cef76d160f29943f9f38a8930e`.

## Fixed-catalog signeR comparison

We refit the release default on all 96 profiles from the earlier signeR panel.
The script verifies the archived input and comparator checksums before fitting.
The comparison uses archived estimates from both external tools.
The profiles span eight cancer exposure distributions, four conditions, and three burdens.
Each clean-profile row represents eight simulations.

| Profiles | MutaGene default | SigProfilerAssignment | signeR | signeR with background |
|---|---:|---:|---:|---:|
| All 96 | 0.19347 | 0.28403 | 0.40704 | 0.33329 |
| Clean, 100 mutations | 0.17152 | 0.24128 | 0.76524 | 0.78840 |
| Clean, 2,000 mutations | 0.09354 | 0.12043 | 0.08803 | 0.08724 |
| Clean, 50,000 mutations | 0.00676 | 0.01026 | 0.00827 | 0.00873 |

All entries report mean total attribution error.
The preceding joint-BIC search scores 0.44720 on clean profiles at 100 mutations.
The current default reverses that earlier ordering against SigProfilerAssignment on the same inputs.
The high-burden refits reproduce the archived joint-BIC errors within numerical tolerance.
The [individual scores](signer-panel-errors.tsv) also include all archived MutaGene baselines.

signeR uses its fixed-catalog mode in both configurations.
Nineteen of twenty evaluable independent-chain pairs exceed 1.01 for maximum active-component R-hat.
These diagnostics limit conclusions about converged signeR performance.
See the [sampling assessment](../signature-benchmark-2026/signer-results/findings.md).

To repeat the refit, restore the files listed in `analysis/local-artifacts.json` that the refit script reads.
Use the environment described below.
Pass the restored repository root and a new output filename.

```sh
/tmp/mutagene-comparison-env/bin/python analysis/release-comparison/refit_signer_panel.py /path/to/restored/repository refit.tsv
```

## Reproduce the paired-signature pilot

Use Python 3.12 and a fresh output directory.
From the repository root, run the following commands.

```sh
uv venv --python 3.12 /tmp/mutagene-comparison-env
uv pip install --python /tmp/mutagene-comparison-env/bin/python . -r analysis/signature-benchmark-2026/requirements-benchmark.txt
/tmp/mutagene-comparison-env/bin/python analysis/release-comparison/compare.py /tmp/mutagene-comparison-run
```

The script writes the reference, counts, SigProfilerAssignment output, and per-profile errors.
The [per-profile scores](scores.tsv) record the generating pair and both errors for each sample.
The [earlier comparator](../signature-benchmark-2026/background-report.qmd) uses empirical exposures and different generation and fitting catalogs.
