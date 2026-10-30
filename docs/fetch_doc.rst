Download reference data
=======================

``mutagene fetch`` downloads example files, cohort data, and reference genome assemblies.
Use ``--help`` with any fetch command to inspect its arguments.

Download example files
----------------------

The examples use ``sample1.maf`` and ``sample2.vcf`` from the MutaGene server.
The command saves both files in the current directory:

.. code-block:: bash

   mutagene fetch examples

Download cohort data
--------------------

Cohort data supplies mutational profiles, observed mutation counts, and sample counts for ranking.
The default source is ``COSMIC``.

.. code-block:: text

   mutagene fetch cohorts [{COSMIC,GDC,MSKCC,ICGC}] [--list | --cohort COHORT]

``--list``, ``-l``
    List cohorts from the selected source.

``--cohort COHORT``, ``-c COHORT``
    Download one cohort by name. This option applies to ``MSKCC``.

List COSMIC cohorts:

.. code-block:: bash

   mutagene fetch cohorts --list

The command reads ``cohorts.tar.gz`` from the current directory and downloads it if necessary.
Use the displayed cohort names with ``mutagene rank --cohort``.

List cohorts from cBioPortal:

.. code-block:: bash

   mutagene fetch cohorts MSKCC --list

Download a reference genome
---------------------------

Select the assembly that matches the coordinates in your mutation file.
MutaGene downloads UCSC reference sequences in 2bit format to ``~/.mutagene/genomes/``.
Repeat an interrupted download command to resume it.

List available assemblies and local files:

.. code-block:: bash

   mutagene fetch genome --list

Download hg19 or hg38:

.. code-block:: bash

   mutagene fetch genome -g hg19
   mutagene fetch genome -g hg38

``--genome``, ``-g``
    Select an assembly using UCSC names, such as ``hg38``, ``hg19``, ``mm10`` or ``mm9``.

``--list``, ``-l``
    List available assemblies and show which files already exist locally.

Pass the full genome path to analysis commands, for example ``-g ~/.mutagene/genomes/hg19.2bit``.
See the `UCSC sequence downloads <http://hgdownload.soe.ucsc.edu/downloads.html>`_ for reference data.
