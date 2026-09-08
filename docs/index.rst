MutaGene Documentation
======================

MutaGene is a Python package for analyzing mutations and mutational processes in cancer.
It provides command-line tools that complement the
`MutaGene website <https://www.ncbi.nlm.nih.gov/research/mutagene/>`_.

How the subcommands fit together
--------------------------------

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

The genome assembly must match the coordinates in the input file; a mismatch is
reported rather than left to produce quiet nonsense. See :doc:`profile_doc`.

Subcommands
-----------

.. toctree::
   :maxdepth: 2

   fetch_doc
   profile_doc
   rank_doc
   motif_doc
   signature_doc
   serve_doc
   common_options

Installation
------------

Requires Python 3.10 or higher::

    pip install mutagene

For the local web interface::

    pip install mutagene[web]

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


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
