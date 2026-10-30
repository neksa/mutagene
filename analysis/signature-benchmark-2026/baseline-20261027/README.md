[![](https://img.shields.io/pypi/v/mutagene.svg)](https://pypi.python.org/pypi/mutagene)
[![](https://readthedocs.org/projects/mutagene/badge/?version=latest)](https://mutagene.readthedocs.io/en/latest/?badge=latest)

# MutaGene

MutaGene is a Python package for analyzing mutations and mutational processes in cancer. It provides command-line tools that complement the [MutaGene website](https://www.ncbi.nlm.nih.gov/research/mutagene/).

## Installation

Requires Python 3.10+.

The source checkout includes features that PyPI 0.9.2.1 does not provide.
Until the next package release, install this checkout with `pip install ".[web]"`.

```bash
pip install mutagene
```

For the local web interface:
```bash
pip install "mutagene[web]"
```

## Subcommands

| Command | Description |
|---------|-------------|
| `mutagene fetch` | Download genomes, cancer datasets, and cohorts from remote sources |
| `mutagene profile` | Create mutational profiles from samples (MAF/VCF input) |
| `mutagene rank` | Rank mutations by driver status using expected mutability |
| `mutagene motif` | Test samples for presence of mutational motifs |
| `mutagene signature` | Decompose profiles into known mutational signatures |
| `mutagene serve` | Start local web server for interactive analysis |

Use `mutagene <command> --help` for detailed usage of each subcommand.

## Quick Start

Try the bundled synthetic example without network access:

```bash
mutagene demo --output-dir mutagene-demo
```

The expected result contains 3 samples, 24 mutations, and 96 channels.
This small example demonstrates execution, not biological inference.

For real mutation files, use a matching reference assembly:

```bash
# Download a genome assembly and the precalculated cohorts
mutagene fetch genome -g hg19
mutagene fetch cohorts
mutagene fetch examples

# Create a mutational profile
mutagene profile -i sample1.maf -g hg19 -o profile.tsv

# List the available cohorts, then rank driver mutations against one
mutagene rank -c
mutagene rank -i sample1.maf -g hg19 --cohort pancancer

# Search for mutational motifs
mutagene motif -i sample1.maf -g hg19 --motif 'C[C>T]G'

# Decompose into COSMIC signatures
mutagene signature -i sample1.maf -g hg19 -s COSMICv2

# Start the web interface
mutagene serve
```

Every subcommand accepts `--params-out FILE` to record the parameters of a run,
and `--params-in FILE` to replay it:

```bash
mutagene profile -i sample1.maf -g hg19 -o profile.tsv --params-out run.json
mutagene profile --params-in run.json -o rerun.tsv
```

The genome assembly must match the coordinates in your input file. Run with `-v`
to see how many reference alleles disagree with the assembly; a large share
means the wrong one was chosen.

## Pipeline and R interchange

Export one sample per column, then fit the labelled SBS96 count matrix:

```bash
mutagene profile -i sample1.maf -g hg19 --per-sample -o matrix.tsv
mutagene signature -i matrix.tsv -f MATRIX -o signatures.tsv --diagnostics-out diagnostics.json
```

Matrix fitting needs no reference genome download.
The first column must contain `Type` and all 96 channel labels.
Remaining columns contain nonnegative integer counts for uniquely named samples.

VCF input supports at most one sample column.
Select the tumor sample and retain its variant calls before analysis.
MutaGene rejects multisample VCFs instead of pooling their records.

## Citation

If you use MutaGene, please cite:

> Goncearenco A, Rager SL, Li M, Sang Q, Rogozin IB, Panchenko AR.
> Exploring background mutational processes to decipher cancer genetic heterogeneity.
> *Nucleic Acids Res.* 2017; 45(W1):W514-W522.
> https://doi.org/10.1093/nar/gkx367

For the driver ranking method (`mutagene rank`):

> Brown AL, Li M, Goncearenco A, Panchenko AR.
> Finding driver mutations in cancer: Elucidating the role of background mutational processes.
> *PLOS Computational Biology* 2019; 15(4): e1006981.
> https://doi.org/10.1371/journal.pcbi.1006981

### Optional annotations and workflow integrations

Fetch an indexed AlphaMissense table with `mutagene fetch alphamissense -g hg38`.
Install the `annotations` extra first.
The [annotation guide](docs/alphamissense.md) explains mapping requirements and the planned driver-model evaluation.
AlphaMissense annotations preserve existing B-scores and driver labels.

The [pipeline examples](examples/pipelines/README.md) cover DVC, Nextflow, and Snakemake.
The [R bridge](examples/r/README.md) exchanges labelled SBS96 matrices through the CLI.
The [controlled benchmark](analysis/signature-benchmark-2026/report.qmd) compares null handling and external attribution.
An [optional SigProfilerAssignment adapter](examples/signatures/README.md) preserves explicit channel alignment.
The [application-note draft](manuscript/README.md) contains the LaTeX source and review requirements.

`--diagnostics-out` reports optimizer status, fallback behavior, and null-component weights for signature fitting.
Python, CLI, and web interfaces now share the repaired `MLE` default.
`MLE-scaled` and `Compat` use the same certified likelihood solver.
The [likelihood guide](docs/likelihood-fitting.md) documents validation, convergence, and compatibility changes.
The [repair reassessment](analysis/signature-benchmark-2026/repaired-report.qmd) compares paired and fresh simulation panels.

Experimental `sparse-bic` and `sparse-aicc` methods search signature subsets with likelihood refits.
The [sparse-fitting guide](docs/sparse-signatures.md) explains the separate treatment of context-independent background.
The [extended-context logo proposal](docs/extended-context-logos.md) describes a separate future feature.
