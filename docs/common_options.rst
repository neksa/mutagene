=============================================
Arguments shared by the analysis subcommands
=============================================

``profile``, ``rank``, ``motif`` and ``signature`` all read variants and all
accept the arguments below.

--------------------------------
1. Variants the caller rejected
--------------------------------

Variant callers record their verdict in a ``FILTER`` field: ``PASS`` for a
variant they stand behind, ``.`` when no filters were applied, and otherwise the
names of the filters it failed, such as ``germline``, ``weak_evidence`` or
``clustered_events``.

**Only PASS variants are counted by default.** A file with no ``FILTER`` column
has rejected nothing, so all of it is used. Pass ``--keep-filtered`` to include
rejected variants.

This matters more than it sounds. A Mutect2 file in which most rows are marked
``germline`` yields a profile that is mostly germline variation rather than
somatic mutation, and the signature decomposition will report that faithfully --
typically as SBS5 with a germline-contamination signature alongside it.

The field is read from the seventh column of a VCF, and from a column named
``FILTER`` in a MAF. GDC's separate ``GDC_FILTER`` column is **not** used: it
mixes quality flags with annotations such as ``NonExonic``, and a non-exonic
mutation is perfectly good input to a mutational profile.

Some MAF converters do not produce a ``FILTER`` column at all and leave the
caller's verdict in a passthrough column whose name carries no meaning. ANNOVAR,
for example, writes Mutect2's verdicts into ``Otherinfo10``. Name that column to
have it honoured::

    mutagene profile -i sample.maf -g hg19 --filter-column Otherinfo10

Such a column is never guessed at: without the option the file is read in full,
because a column called ``Otherinfo10`` announces nothing about its contents.

A column named with ``--filter-column`` and not found in the file is an error
rather than a warning. Carrying on would filter nothing at all, which is the
outcome the option exists to prevent.

-----------------------------
2. Recording a run
-----------------------------

--params-out FILE
    Write the parameters of this run to FILE in JSON format, including the
    resolved genome path, so the run can be reproduced from a different working
    directory.

--params-in FILE
    Take parameters from a FILE written by ``--params-out``. Arguments given on
    the command line take precedence, so a recorded run can be repeated with one
    value changed::

        mutagene profile -i sample.maf -g hg19 -o profile.tsv --params-out run.json
        mutagene profile --params-in run.json -g hg38 -o hg38.tsv

    A file recorded for one subcommand is refused by another.

--------------------------------
3. Exit status
--------------------------------

A run that reads no mutations reports failure rather than writing an empty
output file and exiting 0, so a shell pipeline stops instead of carrying on with
nothing. Finding no motifs, or no signature above threshold, in mutations that
did load is a result and succeeds.
