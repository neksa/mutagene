# Compare signature attribution methods

The controlled benchmark supports SigProfilerAssignment as an optional attribution method.
It does not support a universal replacement of MutaGene's null-signature model.
See the [benchmark report](../../analysis/signature-benchmark-2026/report.qmd) for conditions and limitations.

In a separate environment, install `SigProfilerAssignment==1.1.5`.
Provide a labelled SBS96 count matrix and an explicit signature reference:

```sh
python examples/signatures/sigprofiler_assignment.py \
  --matrix matrix.tsv --reference reference.tsv --output assignment
```

The adapter aligns both tables to SigProfilerAssignment's canonical channel order.
This step matters because the library reorders samples but preserves custom reference row order.
Record reference hashes, package versions, and genome opportunities with the results.
The external method does not reproduce MutaGene's six null components.
Do not compare reconstruction similarity alone when deciding which attribution to trust.
