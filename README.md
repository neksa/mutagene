[![](https://img.shields.io/pypi/v/mutagene.svg)](https://pypi.python.org/pypi/mutagene)
[![](https://img.shields.io/circleci/build/gh/neksa/mutagene?token=e0e77379f7c1c1b136bf15b30494d0a18957e751)](https://circleci.com/neksa/mutagene)
[![](https://readthedocs.org/projects/mutagene/badge/?version=latest)](https://mutagene.readthedocs.io/en/latest/?badge=latest)

------------------
What is MutaGene?
------------------

MutaGene is a set of methods for the analysis of mutations and mutational processes in cancer. The MutaGene Python package
consists of command-line tools that provide direct access to some of the functions available on MutaGene's website https://www.ncbi.nlm.nih.gov/research/mutagene/.

The MutaGene software package includes 5 subpackages: fetch, profile, rank, motif, and signature. Each subpackage has a
unique functionality.


-----------------------------------
How can I install and use MutaGene?
-----------------------------------

**Installation:** MutaGene package is accessible via standard python repository PyPi. The package requires Python3.7 or higher.

It is recommended to use `venv` or `conda` to create an isolated environment for mutagene package: 
```
python3 -m venv env_mutagene
source env_mutagene/bin/activate
```

Then install from pypi
```
pip install mutagene
```

**Usage:** The MutaGene package can be ran on the command line.

```
usage: mutagene [-v] [-V] [-h] {fetch, profile, rank, motif, signature} ...

MutaGene version 0.9.1.0 - Analysis of mutational processes and driver mutations

Global optional arguments:
  -v, --verbose         Print additional messages (-v, -vv)
  -V, --version         Show version and exit
  -h, --help            Show this help message and exit

Choose MutaGene subpackage:
          fetch - Load data such as genomes and cancer datasets from demote sources (alias: download)
          profile - Create a mutational profile given a sample with mutations
          rank - Predict driver mutations by ranking observed mutations with respect to their expected mutability
          motif - Test samples for presence of mutational motifs
          signature - Identify activity of existing mutational signatures in samples or derive new signatures (aliases: identify, decompose)

  {fetch, profile, rank, motif, signature}
```


--------------------------------------
What are MutaGene's five subpackages?
--------------------------------------

1. MutaGene Fetch allows you to download sample data, cohorts, and human genome reference sequences.

2. MutaGene Profile allows you to analyze any set of mutations from one or several cancer samples to identify cancers of
   unknown primary tumor site, to detect the most likely mutational process and to distinguish tumorigenic from normal or benign mutation sets

3. MutaGene Rank allows you to rank mutations with respect to their driver status in a given sample or cohort in a batch mode using pre-calculated or user-provided mutational profiles

4. MutaGene Motif allows you to analyze the presence of mutational motifs in genomic data.

5. MutaGene Signature allows you to analyze the presence of mutational signatures in genomic data.


# module fetch

```
usage: mutagene fetch [-h] {examples,cohorts,genome} ...

Download data from remote repositories and API

optional arguments:
  -h, --help            show this help message and exit

subcommands:
  Choose data source

  {examples,cohorts,genome}
                        additional help available for subcommands
    cohorts             cohorts help
```

# module profile

```
 mutagene profile --help
usage: mutagene profile [-h] [--infile [INFILE ...]] [--outfile [OUTFILE]] [--genome GENOME] [--input-format INPUT_FORMAT]

positional arguments:


optional arguments:
  -h, --help            show this help message and exit
  --infile [INFILE ...], -i [INFILE ...]
                        Input file format
  --outfile [OUTFILE], -o [OUTFILE]
                        Name of output file, will be generated in TSV format
  --genome GENOME, -g GENOME
                        Location of genome assembly file
  --input-format INPUT_FORMAT, -f INPUT_FORMAT
                        Input format: auto, MAF, VCF
```

# module rank

```
usage: mutagene rank [-h] [--infile INFILE] [--genome GENOME] [--outfile [OUTFILE]] [--cohorts-file [COHORTS_FILE]] [--cohort [COHORT]] [--profile PROFILE] [--nsamples NSAMPLES]
                     [--threshold-driver THRESHOLD_DRIVER] [--threshold-passenger THRESHOLD_PASSENGER]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  --infile INFILE, -i INFILE
                        Input file in MAF format
  --genome GENOME, -g GENOME
                        Location of genome assembly file in 2bit format

Optional arguments:
  --outfile [OUTFILE], -o [OUTFILE]
  --cohorts-file [COHORTS_FILE]
                        Location of tar.gz container or directory for cohorts
  --cohort [COHORT], -c [COHORT]
                        Name of cohort with observed mutations

Advanced arguments:
  --profile PROFILE, -p PROFILE
                        Override profile to calculate mutability, may also describe cohort size
  --nsamples NSAMPLES, -n NSAMPLES
                        Override cohort size
  --threshold-driver THRESHOLD_DRIVER, -td THRESHOLD_DRIVER
                        BScore threshold between Driver and Pontential Driver mutations
  --threshold-passenger THRESHOLD_PASSENGER, -tp THRESHOLD_PASSENGER
                        BScore threshold between Pontential Driver and Passenger mutations
```

# module motif
```
usage: mutagene motif [-h] [--infile INFILE] [--genome GENOME] [--input-format {MAF,VCF}] [--motif MOTIF] [--outfile [OUTFILE]] [--window-size WINDOW_SIZE] [--strand {T,N,A,TNA}]
                      [--threshold THRESHOLD] [--save-motif-matches SAVE_MOTIF_MATCHES] [--test {Fisher,Chi2}]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  --infile INFILE, -i INFILE
                        Input file in MAF or VCF format with one or multiple samples
  --genome GENOME, -g GENOME
                        Location of genome assembly file in 2bit format

Optional arguments:
  --input-format {MAF,VCF}, -f {MAF,VCF}
                        Input format: MAF, VCF
  --motif MOTIF, -m MOTIF
                        Motif to search for, use the 'R[C>T]GY' syntax for the motif. Use quotes
  --outfile [OUTFILE], -o [OUTFILE]
                        Name of output file, will be generated in TSV format


Advanced arguments:
  --window-size WINDOW_SIZE, -w WINDOW_SIZE
                        Context window size for motif search, default setting is 50
  --strand {T,N,A,TNA}, -s {T,N,A,TNA}
                        Transcribed strand (T), non-transcribed (N), any (A), or all (TNA default)
  --threshold THRESHOLD, -t THRESHOLD
                        Significance threshold for qvalues, default value=0.05
  --save-motif-matches SAVE_MOTIF_MATCHES
                        Save mutations in matching motifs to a BED file
  --test {Fisher,Chi2}  Statistical test to use

Examples:
# search in sample2.vcf for all preidentified motifs in mutagene using hg19
mutagene motif --infile sample2.vcf --input-format VCF --genome hg19

# search for the presence of the C[A>T] motif in sample1.maf using hg19 not checking for strand-specificity
mutagene motif --infile sample1.maf --input-format MAF --genome hg19 --motif 'C[A>T]' --strand A
```

# module signature
```
usage: mutagene signature [-h] [--infile INFILE] [--genome GENOME] [--signatures {5,10,30,49,53,MGA,MGB,COSMICv2,COSMICv3,KUCAB}] [--input-format {MAF,VCF,TCGI}] [--outfile [OUTFILE]]
                          [--method [METHOD]] [--no-unexplained-variance] [--mutations-threshold MUTATIONS_THRESHOLD] [--keep-only KEEP_ONLY] [--bootstrap]
                          [--bootstrap-replicates BOOTSTRAP_REPLICATES] [--bootstrap-confidence-level BOOTSTRAP_CONFIDENCE_LEVEL] [--bootstrap-method {t,p}]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  --infile INFILE, -i INFILE
                        Input file in VCF or MAF format
  --genome GENOME, -g GENOME
                        Location of genome assembly file in 2bit format
  --signatures {5,10,30,49,53,MGA,MGB,COSMICv2,COSMICv3,KUCAB}, -s {5,10,30,49,53,MGA,MGB,COSMICv2,COSMICv3,KUCAB}
                        Collection of signatures to use

Optional arguments:
  --input-format {MAF,VCF,TCGI}, -f {MAF,VCF,TCGI}
                        Input format: MAF, VCF
  --outfile [OUTFILE], -o [OUTFILE]
                        Name of output file, will be generated in TSV format

Advanced arguments:
  --method [METHOD], -m [METHOD]
                        Method defines the function minimized in the optimization procedure
  --no-unexplained-variance, -U
                        Do not account for unexplained variance (non-context dependent mutational processes and unknown signatures)
  --mutations-threshold MUTATIONS_THRESHOLD, -t MUTATIONS_THRESHOLD
                        Only report signatures with mutations above the threshold
  --keep-only KEEP_ONLY, -k KEEP_ONLY
                        Keep only the signatures in the list, separated by commas e.g. 1,3,5

Bootstrap-specific arguments:
  --bootstrap, -b       Use the bootstrap to calculate confidence intervals
  --bootstrap-replicates BOOTSTRAP_REPLICATES, -br BOOTSTRAP_REPLICATES
                        Number of bootstrap replicates
  --bootstrap-confidence-level BOOTSTRAP_CONFIDENCE_LEVEL, -bcl BOOTSTRAP_CONFIDENCE_LEVEL
                        Confidence level
  --bootstrap-method {t,p}, -bm {t,p}
                        Bootstrap method (t: t-distribution, p: percentile)
```

------------------------
How can I cite MutaGene?
------------------------
If you use MutaGene, you should cite:

Goncearenco A, Rager SL, Li M, Sang Q, Rogozin IB, Panchenko AR
Exploring background mutational processes to decipher cancer genetic heterogeneity. Nucleic Acids Res. 2017; 45(W1):W514–W522.
https://doi.org/10.1093/nar/gkx367

Additionally, if you use the MutaGene rank subpackage, you should cite the driver prediction method as:

Brown AL, Li M, Goncearenco A, Panchenko AR Finding driver mutations in cancer: Elucidating the role of background mutational processes.
PLOS Computational Biology 2019; 15(4): e1006981. https://doi.org/10.1371/journal.pcbi.1006981
