[![](https://img.shields.io/pypi/v/mutagene.svg)](https://pypi.python.org/pypi/mutagene)
[![](https://readthedocs.org/projects/mutagene/badge/?version=latest)](https://mutagene.readthedocs.io/en/latest/?badge=latest)

# MutaGene

MutaGene is a Python package for analyzing mutations and mutational processes in cancer. It provides command-line tools that complement the [MutaGene website](https://www.ncbi.nlm.nih.gov/research/mutagene/).

## Installation

Requires Python 3.10+.

```bash
pip install mutagene
```

For the local web interface:
```bash
pip install mutagene[web]
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

```bash
# Download a genome assembly and the precalculated cohorts
mutagene fetch genome hg19
mutagene fetch cohorts

# Create a mutational profile
mutagene profile -i sample.maf -g hg19 -o profile.tsv

# List the available cohorts, then rank driver mutations against one
mutagene rank -c
mutagene rank -i sample.maf -g hg19 --cohort Lung_Adenocarcinoma

# Search for mutational motifs
mutagene motif -i sample.maf -g hg19 --motif 'C[C>T]G'

# Decompose into COSMIC signatures
mutagene signature -i sample.maf -g hg19 -s COSMICv2

# Start the web interface
mutagene serve
```

Every subcommand accepts `--params-out FILE` to record the parameters of a run,
and `--params-in FILE` to replay it:

```bash
mutagene profile -i sample.maf -g hg19 -o profile.tsv --params-out run.json
mutagene profile --params-in run.json -o rerun.tsv
```

The genome assembly must match the coordinates in your input file. Run with `-v`
to see how many reference alleles disagree with the assembly; a large share
means the wrong one was chosen.

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
