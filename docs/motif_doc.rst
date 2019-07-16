==============================
Description
==============================

----------------------------------------
Functionality
----------------------------------------
Use mutagene motif to search for the presence of mutational motifs in mutational data

----------------------------------------
Motif Description
----------------------------------------

A motif is a characteristic pattern of DNA change, often associated with a carcinogen or a biological process. 
A motif's pattern includes a characteristic DNA context and a single-base substitution. 

----------------------------------------
Motif Representation in MutaGene
----------------------------------------

MutaGene represents motifs in a string of characters, where the characters not in brackets represent the unmutated DNA context,
and the characters in brackets represent the single-base substitution. The motif must be in quotes to be recognized by MutaGene.

Motif Examples
--------------

"A[C>A]G" represents the DNA sequence "ACG" mutated into the DNA sequence "AAG"

"C[G>T]" represents the DNA sequence "CG" mutated into the DNA sequence "CT"

"[C>A]C" represents the DNA sequence "CC" mutated into the DNA sequence "AC"

--------------------------------------------------
Further Reading on Motifs
--------------------------------------------------

The publication `Mutational signatures and mutable motifs in cancer genomes <https://doi.org/10.1093/bib/bbx049>`_ describes motifs and their uses

==============================
Motif Subpackage Documentation
==============================

Option 1 - Search for motifs in genomic data

Command: mutagene motif search [--infile][--genome][--motif][--outfile][--window-size][--strand][--help]

Required arguments (must be specified):

--infile, -i
     Input file in MAF or VCF format with one or multiple samples

--genome, -g
    Location of genome assembly file in 2bit format

Optional arguments (can be specified):

--motif, -m
    Motif to search for, use the 'R[C>T]GY' syntax for the motif. Use quotes. If left unspecified, MutaGene will search
    for all pre-identified motifs

--outfile, -o
    Name of output file, will be generated in TSV format

--window-size, -w
    Context window size for motif search, default setting is 50

--strand, -s
    Transcribed strand (+), non-transcribed (-), any (*default), or all (-+*)

--help, -h
    show this help message and exit

Examples:

Option 1 - Search for motifs in mutational data

To search for presence of the C[A>T] motif in sample1.maf using hg19

    $ mutagene motif search -i sample1.maf -g hg19 -m 'C[A>T]'

To search sample2.vcf for all preidentified motifs in mutagene using hg18

    $ mutagene motif search -i sample2.vcf -g hg18

To search sample1.maf for all preidentified motifs in mutagene using hg18 and have output go into file named motif_results

    $ mutagene motif search -i sample1.maf -g hg18 -o motif_results

Option 2 - List all pre-identified motifs in MutaGene

    $ mutagene motif list

=============================
How to Interpret Motif Output
=============================

If no motifs are significantly present in the data, the output will say: "WARNING No significant motif matches found"

If the presence of a motif is significant in the data, the output will show a table with the following headers:

======  ======  =========   ===========  ================  ===========  ===================  ===================
Sample   Name     Motif       Strand       Enrichment        pvalue      mutations_low_est    mutations_high_est
======  ======  =========   ===========  ================  ===========  ===================  ===================

- Sample: Name of Sample. If input file contains multiple samples, output will be stratified per sample.

- Name: Name of motif. If -m/--motif argument is given, name will be "Custom motif".

- Motif: Motif searched for in data

- Strand: DNA Strand that motif was searched for on. '+': transcribed strand, '-': non-transcribed strand, "*": any strand,
  or "-+*": all strands

- Enrichment: Quantitative measure of motif's prevalence, significant if greater than one.

- pvalue: Fisher's p-value for motif significance

- mutations_low_est: Conservative estimate for number of mutations (of total number in input file) that match motif

- mutations_high_est: Maximum number of mutations (of total number in input file) that match the motif

Sample Output and Interpretation (from running mutagene motif search -i sample1.maf -g hg19):

+------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
| sample                       | name       | motif      | strand | enrichment        | pvalue                 | mutations_low_est | mutations_high_est |
+------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
| TCGA-50-6593-01A-11D-1753-08 | C>T in CpG | [C>T]G     | '*'    | 4.586718025481874 | 1.0181609110804669e-06 | 15                | 18.0               |
+------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+

In sample "TCGA-50-6593-01A-11D-1753-08", an estimated 15-18 mutations are thought to be contributed by the
mutagenic process(es) represented by the C>T in CpG motif ([C>T]G). The measures of significance, enrichment and Fisher's pvalue calculations,
are 4.586718025481874 and 1.0181609110804669e-06 respectively.