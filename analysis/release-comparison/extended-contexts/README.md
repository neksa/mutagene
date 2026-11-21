# Fit pentanucleotide and indel signatures

The [figure](contexts.pdf) shows the leading reference signature in each example.
The [table](example-fits.tsv) gives exposures and estimated mutation counts from the default sparse BIC fitter.
The script chooses the first sample in each source matrix before fitting.

| Context | Sample | Observations | Catalog | Highest fitted exposure |
|---|---|---:|---|---|
| SBS1536 | PD10010a | 2,430 substitutions | 60 PCAWG SignatureAnalyzer profiles | SBS39_P, 14.0% |
| ID83 | LUAD-US_SP50263 | 522 indels | 23 COSMIC v3.4 GRCh37 profiles | ID1, 34.3% |

The SBS1536 fit retains both flanking bases on either side of each substitution.
Each of the 96 inner channels splits into 16 outer-base combinations.
Panel A shows all 1,536 probabilities. Panel B expands the inner channel with the largest total probability.
The axes identify the 5′ and 3′ outer bases. Cell values give percent of total signature mass.
Panel C shows the five most probable ID83 channels in ID1.
The class labels encode length, event type, repeat or microhomology class, and bin.
Fitted exposures describe samples. Reference profiles describe the mutation distribution within a signature.
The example fits have no independent exposure ground truth.

The [runner](example.py) downloads all four inputs from pinned commits.
It aligns channel labels and runs the public CLI for each sample.
It converts the SignatureAnalyzer annotation columns to a labelled reference matrix.
It saves the prepared references, sample matrices, diagnostics, and exposure tables in the work directory.
The runner uses the [SignatureAnalyzer SBS reference](https://github.com/getzlab/SignatureAnalyzer/blob/f205be05f6e180bd04c3beaffc354b73ae3ccf9f/signatureanalyzer/ref/PCAWG/sa_PCAWG_sbs.tsv)
and [breast cancer counts](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator/blob/7bbd3dd5dfff7d90de28e091cc2e2c98472f1512/SigProfilerMatrixGenerator/references/matrix/BRCA_bench/SBS/BRCA_bench.SBS1536.all).
The indel example uses the [COSMIC v3.4 GRCh37 catalog](https://github.com/SigProfilerSuite/SigProfilerAssignment/blob/308ea01dce1ba1136aac6f2246f535284b5234d8/SigProfilerAssignment/data/Reference_Signatures/GRCh37/COSMIC_v3.4_ID_GRCh37.txt)
and [protocol ID83 matrix](https://github.com/SigProfilerSuite/Protocol_SigProfiler/blob/25e7b10bb65b871363a25ebda69fe7101c91e578/input_example/matrices/LUAD.ID83.all).

Use the benchmark Python environment from the repository root.

```sh
PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python analysis/release-comparison/extended-contexts/example.py /tmp/mutagene-contexts
```

## Known-mixture recovery

The [pilot](recovery.py) fits 60 synthetic profiles through the CLI.
For each context, the first six published reference columns define the fitting catalog.
Each of the 15 pairs generates one 70:30 mixture per burden. The earlier column supplies 70%.
The sampler uses seed 7321 and multinomial draws.
The score is half the sum of absolute exposure errors across all six columns.

| Mean attribution error | 100 mutations | 1,000 mutations |
|---|---:|---:|
| SBS1536 | 0.0598 | 0.0158 |
| ID83 | 0.0421 | 0.0129 |

Each mean summarizes 15 profiles. The [individual scores](recovery.tsv) retain the generating pair and mutation burden.
This pilot tests clean mixtures within a small known catalog.
It measures neither missing-signature robustness nor performance against other tools.
After preparing the inputs with the example runner, execute the following command from the repository root.

```sh
PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python analysis/release-comparison/extended-contexts/recovery.py /tmp/mutagene-contexts
```
