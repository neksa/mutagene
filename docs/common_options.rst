Shared analysis options
=======================

``profile``, ``rank``, ``motif`` and ``signature`` accept these options when reading variants.

Filter variants
---------------

MutaGene accepts ``PASS``, ``.`` and empty values in the ``FILTER`` field.
It excludes other values unless you specify ``--keep-filtered``.
If the file has no ``FILTER`` column, MutaGene reads all variants.

Variant callers use ``PASS`` for accepted variants and ``.`` when they apply no filters.
Other values name failed filters, such as ``germline``, ``weak_evidence`` or ``clustered_events``.

.. mermaid::

   flowchart TD
       V["variant row"] --> K{"--keep-filtered<br/>given?"}
       K -- yes --> USE["counted"]
       K -- no --> C{"FILTER column<br/>present?"}
       C -- "no column" --> USE
       C -- yes --> F{"value is<br/>PASS, . or empty?"}
       F -- yes --> USE
       F -- no --> DROP["excluded,<br/>and reported"]

       classDef keep fill:#e8f4ea,stroke:#4a7c59,color:#1b3a29
       classDef drop fill:#fdeaea,stroke:#b1483f,color:#4a1a16
       class USE keep
       class DROP drop

MutaGene reads the seventh VCF column or the MAF column named ``FILTER``.
It does not use ``GDC_FILTER``, which combines quality flags with annotations such as ``NonExonic``.
Non-exonic mutations remain valid input for mutational profiles.

Some converters preserve filter values under another column name.
For example, ANNOVAR stores Mutect2 filter values in ``Otherinfo10``.
Specify that column explicitly:

.. code-block:: bash

   mutagene profile -i sample.maf -g hg19 --filter-column Otherinfo10

MutaGene does not infer the meaning of arbitrary column names.
If the requested column does not exist, the command reports an error.

Record and repeat a run
-----------------------

``--params-out FILE``
    Save parameters as JSON, including the resolved genome path.
    This path lets you repeat the run from another working directory.

``--params-in FILE``
    Read saved parameters. Explicit command-line arguments take precedence.
    Each subcommand accepts only its own parameter files.

Save a run, then repeat it with another assembly and output path:

.. code-block:: bash

   mutagene profile -i sample.maf -g hg19 -o profile.tsv --params-out run.json
   mutagene profile --params-in run.json -g hg38 -o hg38.tsv

Interpret exit status
---------------------

A run fails if it reads no mutations.
Finding no motifs or no signatures above the reporting threshold remains a valid result when mutations are available.
