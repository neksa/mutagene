[![](https://img.shields.io/pypi/v/mutagene.svg)](https://pypi.python.org/pypi/mutagene)
[![](https://img.shields.io/circleci/build/gh/neksa/mutagene?token=e0e77379f7c1c1b136bf15b30494d0a18957e751)](https://circleci.com/neksa/mutagene)
[![](https://readthedocs.org/projects/mutagene/badge/?version=latest)](https://mutagene.readthedocs.io/en/latest/?badge=latest)
[![](https://pyup.io/repos/github/neksa/mutagene/shield.svg)](https://pyup.io/repos/github/neksa/mutagene/)

------------------
What is MutaGene?
------------------

MutaGene is a set of methods for the analysis of mutations and mutational processes in cancer. The MutaGene Python package
consists of command-line tools that provide direct access to some of the functions available on MutaGene's website https://www.ncbi.nlm.nih.gov/research/mutagene/.

The MutaGene software package includes 5 subpackages: fetch, profile, rank, motif, and signature. Each subpackage has a
unique functionality.

--------------------------------------
What are MutaGene's five subpackages?
--------------------------------------

| 1. MutaGene Fetch allows you to download sample data, cohorts, and human genome reference sequences.

| 2. MutaGene Profile allows you to analyze any set of mutations from one or several cancer samples to identify cancers of
   unknown primary tumor site, to detect the most likely mutational process and to distinguish tumorigenic from normal or benign mutation sets

| 3. MutaGene Rank allows you to rank mutations with respect to their driver status in a given sample or cohort in a batch mode using pre-calculated or user-provided mutational profiles

| 4. MutaGene Motif allows you to analyze the presence of mutational motifs in genomic data.

| 5. MutaGene Signature allows you to analyze the presence of mutational signatures in genomic data.

-----------------------------------
How can I install and use MutaGene?
-----------------------------------

**Installation:** MutaGene package is accessible via standard python repository PyPi. The package requires Python3.6 or higher.

**Usage:** The MutaGene package can be ran on the command line.

Here are `detailed instructions <https://www.ncbi.nlm.nih.gov/research/mutagene/package>`_ to install and use MutaGene.
- add links to command-line MutaGene documentation and, of course, subpackages

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
