MutaGene documentation
======================

MutaGene analyzes mutations and mutational processes in cancer.
Use the command line for repeatable analyses or the local web interface to submit files and inspect results.
The package complements the `MutaGene website <https://www.ncbi.nlm.nih.gov/research/mutagene/>`_.

Install MutaGene
----------------

MutaGene requires Python 3.10 or higher.

.. code-block:: bash

   pip install mutagene

For the local web interface, install the web dependencies:

.. code-block:: bash

   pip install mutagene[web]

Choose a task
-------------

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - Task
     - Guide
   * - Download reference data
     - :doc:`fetch_doc`
   * - Count mutations by sequence context
     - :doc:`profile_doc`
   * - Estimate signature exposures
     - :doc:`signature_doc`
   * - Rank potential driver mutations
     - :doc:`rank_doc`
   * - Test mutational motifs
     - :doc:`motif_doc`
   * - Submit files in a browser
     - :doc:`serve_doc`

Match the reference genome to your input coordinates before analysis.
See :doc:`profile_doc` for reference-allele checks.

How the commands connect
------------------------

.. mermaid::

   flowchart LR
       MAF["MAF / VCF<br/>mutations"] --> P["profile"]
       G[("2bit genome<br/>assembly")] -.-> P
       P --> PROF["96-channel<br/>profile"]
       PROF --> S["signature"]
       PROF --> R["rank"]
       MAF --> R
       MAF --> M["motif"]
       G -.-> R
       G -.-> M
       S --> EXP["signature<br/>exposures"]
       R --> DRV["ranked<br/>driver mutations"]
       M --> ENR["motif<br/>enrichment"]

       classDef out fill:#e8f4ea,stroke:#4a7c59,color:#1b3a29
       classDef cmd fill:#e6eefc,stroke:#3b5ea8,color:#16264a
       class P,S,R,M cmd
       class PROF,EXP,DRV,ENR out

.. toctree::
   :hidden:
   :maxdepth: 1

   fetch_doc
   profile_doc
   signature_doc
   rank_doc
   motif_doc
   serve_doc
   common_options

Citation
--------

If you use MutaGene, please cite:

Goncearenco A, Rager SL, Li M, Sang Q, Rogozin IB, Panchenko AR
Exploring background mutational processes to decipher cancer genetic heterogeneity.
*Nucleic Acids Res.* 2017; 45(W1):W514-W522.
https://doi.org/10.1093/nar/gkx367

For the driver ranking method (``mutagene rank``):

Brown AL, Li M, Goncearenco A, Panchenko AR
Finding driver mutations in cancer: Elucidating the role of background mutational processes.
*PLOS Computational Biology* 2019; 15(4): e1006981.
https://doi.org/10.1371/journal.pcbi.1006981


Search and indices
------------------

* :ref:`search`
* :ref:`genindex`
* :ref:`modindex`
