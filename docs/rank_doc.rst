===============
Description
===============

Use mutagene rank to analyze genes and compare observed mutational frequencies to expected background mutability to identify potential drivers.
Cancer-specific mutational models can be used.

To read more about the MutaGene's ranking method, read the publication
`Finding driver mutations in cancer: Elucidating the role of background mutational processes <https://doi.org/10.1371/journal.pcbi.1006981>`_

==================================
Rank Subpackage Documentation
==================================

Command: mutagene rank [-h] [--infile [INFILE [INFILE ...]]] [--outfile [OUTFILE]] [--genome GENOME]
[--cohorts-file COHORTS_FILE] [--cohort COHORT] [--profile PROFILE] [--nsamples NSAMPLES]
[--threshold-driver THRESHOLD_DRIVER] [--threshold-passenger THRESHOLD_PASSENGER]

Required Arguments (must be specified):

--cohort, -c
    Name of cohort with observed mutations

--infile, -i
    Input file in MAF format

--genome, -g
    Location of genome assembly file in 2bit format

Optional Arguments (can be specified):

--outfile, -o
    Name of output file, will be generated in TSV format

--cohorts-file
                    Location of tar.gz container or directory for cohorts

--cohort, -c
                    Name of cohort with observed mutations

--profile, -p
                    Override profile to calculate mutability, may also
                    describe cohort size

--nsamples, -n
                    Override cohort size

--threshold-driver, -td
                    BScore threshold between Driver and Pontential Driver
                    mutations

--threshold-passenger, -tp
                    BScore threshold between Pontential Driver and
                    Passenger mutations
--help, -h
    show this help message and exit

Examples:

To rank mutations in sample1.maf using hg19 and gcb_lymphomas cohort:

    $ mutagene rank -i sample1.maf -g hg19 -c gcb_lymphomas

----------------------------------------------------
Cohorts Available for Analysis
----------------------------------------------------

=================================
How to Interpret Rank Output
=================================

The output will show a table with the following headers:

======  =============  =============  ================  ===========  ===================  =========
 gene    mutation        mutability    observed           bscore      qvalue                label
======  =============  =============  ================  ===========  ===================  =========

- gene:

- mutation:

- mutability: expected mutation rate in a particular DNA context

- observed: observed mutational frequencies

- bscore: a binomial p-value for the observed number of occurences of mutation in comparison to the expected
  mutability that is defined by the local DNA context of the mutated nucleotide.

- qvalue: Bscore corrected for multiple testing with Benjamini-Hochberg FDR method

- label: Prediction of cancer drivers, Potential drivers, and Passengers is based on the thresholds established
  for the Bscore optimized using this benchmark datasets.

Sample Output and Interpretation:
