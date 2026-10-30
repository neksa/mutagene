# Export sample count matrices

The profile command exports one sample per column with `--per-sample`.
Provide a reference genome that matches the mutation coordinates.

```sh
mutagene profile -i cohort.maf.gz -g hg38 --per-sample -o matrix.tsv
```

The first column contains `Type` and the 96 SBS channel labels.
Each remaining column contains nonnegative integer counts for one uniquely named sample.
The reader restores the canonical channel order from labels.
Duplicate sample names across input files cause an error.

The profile and motif commands accept plain text, gzip files, and standard input with `-i -`.
The input helper stores standard input in a seekable temporary stream.

VCF input supports at most one sample column.
Select one tumor sample before analysis.
The readers exclude records whose genotype contains only reference alleles.
