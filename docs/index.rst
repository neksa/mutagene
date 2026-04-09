MutaGene Documentation
======================

MutaGene is a Python package for analyzing mutations and mutational processes in cancer.
It provides command-line tools that complement the
`MutaGene website <https://www.ncbi.nlm.nih.gov/research/mutagene/>`_.

Subcommands
-----------

.. toctree::
   :maxdepth: 2

   fetch_doc
   profile_doc
   rank_doc
   motif_doc
   signature_doc

Installation
------------

Requires Python 3.8 or higher::

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
