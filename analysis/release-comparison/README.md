# Compare the release default with SigProfilerAssignment

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

## Reproduce

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
