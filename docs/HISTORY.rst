Development
===========

Please use the `MutaGene repository <https://github.com/neksa/mutagene/>`_.
Pull requests gladly accepted.
Issues should be reported at the github issue tracker.

Running tests
-------------

Please check the tests by running them with::

    python setup.py test

New features should have test code sent with them.

Changes
=======

0.8.6.6 Bug-Fix Release
-----------------------
* BUGFIX: Mutational profile was not incorrectly calculated for MAF files with multiple samples which affected decomposition for COMIC 30 and 49 signature sets

0.8.6.5 Bug-Fix Release
-----------------------
* MAF file loading improved for GDC and MSKCC data sources. More meaningfull error messages

0.8.6.4 Release
-----------------------
* testing and development releases are not available in pip mirrors, bumping version

0.8.6.4dev1 Testing release
-----------------------

* added handling of VCF files to motif analysis
* Signatures from COSMIC v3 available as signature set "49"


0.8.5.1 Bug-Fix Release
-----------------------

* several error messages downgraded in log level
* correct handling of missing parameters

0.8.5 Release
-------------

* Functionality available 'fetch_genomes', 'fetch_cohorts', and 'rank'
* Uploaded to PyPi as mutagene


Contributions
=============

Project started by Alexander Goncearenco @neksa in 2011 in Anna Panchenko research lab.
Contributions from Caroline Cunningham, Stephanie Rager, and Anna-Leigh Brown.

This project was supported by Intramural Research Grant of the NLM, NIH.
