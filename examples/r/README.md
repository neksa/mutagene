# Exchange SBS96 matrices with R

The bridge sends labelled integer count matrices to the MutaGene CLI.
It preserves sample names and returns an ordinary R data frame.
It requires an installed `mutagene` executable on `PATH`.
It does not require `reticulate` or embed Python in R.

```r
source("examples/r/mutagene_bridge.R")
counts <- as.matrix(read.delim("matrix.tsv", row.names=1, check.names=FALSE))
exposures <- fit_mutagene(counts, "signatures.tsv")
```

`MutationalPatterns::mut_matrix()` produces a compatible SBS96 count matrix.
Supply mutation counts, not normalized fractions or trinucleotide opportunity corrections.
Preserve the original channel labels and sample columns.
Check genome assembly and reference opportunities before comparing WGS and WES estimates.
This bridge does not constitute a Bioconductor package or validate cross-tool scientific equivalence.

See the [MutationalPatterns documentation](https://bioconductor.org/packages/MutationalPatterns/).
