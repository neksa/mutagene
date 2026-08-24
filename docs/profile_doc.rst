==================================
Profile: mutational profiles
==================================

-----------------
1. Description
-----------------

``mutagene profile`` counts the mutations in one or more samples by their
substitution type and immediate sequence context, producing a 96-channel
mutational profile. The profile is the input to signature decomposition and to
the background mutability model used by ``mutagene rank``.

Reference context is read from a 2bit genome assembly, so the assembly must
match the coordinates in the input file. See section 4.

-----------------
2. Arguments
-----------------

**Command:** ``mutagene profile [arguments]``

Required arguments:

--infile, -i
    Input file in MAF or VCF format, with one or multiple samples. May be given
    more than once; the profiles are pooled.

--genome, -g
    Genome assembly: a name such as ``hg19``, or the path to a 2bit file.
    A bare name is resolved against the directory used by ``mutagene fetch``.

Optional arguments:

--outfile, -o
    Name of the output file, in TSV format. Defaults to the screen.

--input-format, -f
    ``auto`` (default), ``MAF`` or ``VCF``. Automatic detection reads the header
    to decide.

--keep-filtered
    Include variants the caller rejected. See section 5.

--params-out FILE
    Write the parameters of this run to FILE as JSON.

--params-in FILE
    Take parameters from a FILE written by ``--params-out``. Arguments given on
    the command line take precedence.

--help, -h
    Show the help message and exit.

------------------------------------
3. How to interpret profile output
------------------------------------

The output has two columns: the mutation channel, and the number of mutations
in the input matching it.

Running ``mutagene profile -i sample2.vcf -g hg19`` produces output beginning::

    A[C>A]A	119
    A[C>G]A	174
    A[C>T]A	475
    A[T>A]A	53

So 119 mutations in ``sample2.vcf`` are C>A substitutions with an A on either
side, 174 are C>G in the same context, and so on.

**Channels are pyrimidine-centric.** Every substitution is expressed with C or T
as the reference base; a G>A mutation is counted as its reverse-complement
C>T, with the flanking bases complemented and swapped. The 96 channels are the
6 substitution types (C>A, C>G, C>T, T>A, T>C, T>G) in each of the 16
5'-3' contexts.

**Row order is not alphabetical.** Rows follow 5' base, then 3' base, then
substitution type -- the order given by ``get_profile_attributes_dict()``, which
is also the row order of the signature matrices that profiles are decomposed
against. Sorting the channel labels alphabetically produces a *different*
permutation of the same 96 channels. Code that pairs a profile with a signature
matrix must not do so by sorted key: the result looks plausible and is wrong.

-----------------------------------
4. Choosing the genome assembly
-----------------------------------

The reference base at each position is looked up in the assembly and compared
with the reference allele in the input. An assembly that does not match the
coordinates produces disagreement at a large fraction of positions, and those
mutations are dropped for want of a usable context.

Run with ``-v`` to see the counts. A high proportion of mismatches is reported
explicitly, for example::

    WARNING 2374 mutations have a reference allele that matches neither strand
    ERROR Found 2374 reference allele mismatches (53% of mutations examined).
          This suggests the wrong genome assembly was selected. Please try hg38.

A reference allele reported on the opposite strand is normal and is handled, not
counted as a mismatch.

------------------------------------
5. Variants the caller rejected
------------------------------------

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

Some MAF converters do not produce a ``FILTER`` column at all, leaving the
caller's verdict in a passthrough column under another name. Nothing can be
honoured in that case; check for such a column before trusting the counts.
