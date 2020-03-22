

# Identify:
    - bootstrap for MAF samples
    - How to make bootstrap results compatible with non-bootstrap results
    - Limit identify to a subset of signatures
    - Create subsets of signatures based on correlation clustering for de-correlation and avoiding collinearity
    - Create named signature sets and named subsets
    - clustergrammer for reporting and visualization of results of identify/decomposition
    - How to output reports?
    - rpy2 or reticulate for integration with R reports
    - BayesianOptimization for global optimization

    - restore functionality of benchmark module

# Motifs
    + output motif matches in a separate file
    - Motif calibration depending on window size, motif length and the number of mutations
    - Comparison of motifs with signatures

# Ranking
    - Implement ranking of nucleotide mutations from VCF
    - Precalculate cohorts
    - rank should accept either MAF (or VCF) or precalculated data as input
    - Profile, Cohort size and observed mutations is required input. By default they are all coming from the same source
    - Define command line use cases
    - Rename cohort to data bundle: samples, 

# Fetch 
    - test MSKCC
    - test GDC
    - test ICGC
    - consider other sources?
    - convert genomes to 2bit
    - Autodownloading genome assembly
    - Re-enable genome-less version (via ENSEMBL) for testing a limited number of mutations or even a single mutation

# General
    - Autodection of genome assembly from VCF and MAF
    - Upgrade website to use mutagene library
    - Accept YAML config file as input for reproducibility
    - Create CWL wrappers and dockerized container versions for pipline deployment
    - Write a parseable human readable log file containing provisioning information for reproducibility. Should include all parameters, versions / config file that can be used as an input

# Tests
    + currently package test suite is failing
    - need to make sure we have a CLI test for every function (not every option combination though)
    ? what other test files can we bundle (licensing issues), perhaps cell lines?

# Documentation
    * readthedocs documentation
        Example: https://plastid.readthedocs.io/en/latest/quickstart.html

    * Improve README
     - main section
     - Installation
     - Changes
     - Documentation
     - Usage (printout of argparse usage)

    * Examples:
        - cell lines examples

# Links:

    - https://github.com/mskcc/vcf2maf
    - https://github.com/PoisonAlien/maftools
    - https://riptutorial.com/bioinformatics/example/14312/mutation-annotation-format--maf-

    - http://docs.h5py.org/en/latest/high/dataset.html#dataset