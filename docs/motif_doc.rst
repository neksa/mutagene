==============================
Description
==============================

----------------------------------------
Functionality
----------------------------------------
Use "mutagene motif" to search for the presence of mutational motifs in mutational data

----------------------------------------
Motif Definition
----------------------------------------

A motif is defined as a characteristic pattern of DNA mutation and its local DNA context. It is often associated with a specific carcinogen or a biological process.

----------------------------------------
Motif Representation
----------------------------------------

MutaGene represents motifs as a string of DNA characters, where characters in brackets represent the single-base substitutions and characters outside brackets represent the unmutated DNA context. The motif must be in quotes to be recognized by MutaGene.

Motif Examples
--------------

"A[C>A]G" represents the DNA sequence "ACG" mutated into the DNA sequence "AAG"

"C[G>T]" represents the DNA sequence "CG" mutated into the DNA sequence "CT"

"[C>A]C" represents the DNA sequence "CC" mutated into the DNA sequence "AC"

--------------------------------------------------
Further Reading on Motifs
--------------------------------------------------

The publication `Mutational signatures and mutable motifs in cancer genomes <https://doi.org/10.1093/bib/bbx049>` describes motifs and their uses

==============================
Motif Subpackage Documentation
==============================

Option 1 - Search for motifs in genomic data

Command: mutagene motif search [--infile][--genome][--motif][--outfile][--window-size][--strand][--help]

Required arguments (must be specified):

--infile, -i
     Input file in MAF format with one or multiple samples

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
    Search for the presence of a motif on Transcribed strand (+), non-transcribed (-), any strand (=) - default option, or all strands (-+=)

--threshold, -t
    Significance threshold for qvalue (corrected p-value) reporting, default value=0.05

--help, -h
    show this help message and exit

------------------------------------
Window Size Parameter Explanation
------------------------------------
MutaGene counts window size as the number of DNA bases searched for from the first base of the DNA sequence gathered up to but not including the mutated base

Therefore, the effective length of the DNA sequence searched is 2 * window-size + 1

It may be advantageous to use a window size longer than the default 50 bases if the motif is longer than three nucleotides,
as this motif is likely to appear less frequently in the DNA context. Similarly, if the motif is shorter than three nucleotides,
it may be advantageous to use a window size shorter than the default 50 bases, as the motif is likely to appear in DNA more frequently.

------------------------------------
Strand Parameter Explanation
------------------------------------
MutaGene can search for the presence of a motif on the transcribed or non-transcribed DNA strands or both strands.
This information is gathered from the input file provided by the user.

Analyzing for the presence on a transcribed or non-transcribed strand is advantageous when a mutational process is
known to have mutations with a transcriptional strand bias. For instance, the APOBEC1/3A/B family is known to be
associated with mutational processes that have a transcriptional strand bias of mutations in exons.

The transcription strand refers to the coding DNA strand, and the non-transcription strand refers to the template DNA strand.

=============================
How to Interpret Motif Output
=============================

If no motifs are significantly present in the data, the output will say: "WARNING No significant motif matches found"

If the presence of a motif is significant in the data, the output will show a table with the following headers:

======  ======  =========   ===========  ================  ===================  ===================  ===========  ===========
Sample   Name     Motif       Strand       Enrichment       mutations_low_est    mutations_high_est   pvalue       qvalue
======  ======  =========   ===========  ================  ===================  ===================  ===========  ===========

- Sample: Name of Sample. If input file contains multiple samples, output will be stratified per sample.

- Name: Name of motif. If -m/--motif argument is given, name will be "Custom motif".

- Motif: Motif searched for in data

- Strand: DNA Strand that motif was searched for on. '+': transcribed strand, '-': non-transcribed strand, "=": any strand, "+-=": all strands.

- Enrichment: Quantitative measure of motif's prevalence, significant if greater than one.

- mutations_low_est: Conservative estimate for number of mutations (of total number in input file) that match motif

- mutations_high_est: Maximum number of mutations (of total number in input file) that match the motif

- pvalue: Fisher's p-value for motif significance

- qvalue: Fisher's p-value with Benjamini-Hochberg correction for motif significance

----------------------------------
How to Interpret Enrichment Output
----------------------------------
Enrichment is modeled off of a risk ratio, meaning that a motif’s enrichment is essentially a ratio between the
probability of a motif appearing in a cancer sample’s DNA mutations and the probability of a motif appearing in a
cancer sample’s DNA context.

Because enrichment is modeled off a risk ratio, it can be interpreted the same way. The result of enrichment minus one
is the percent overrepresentation of a motif. For example, if enrichment is 1.5, it means that there is a 50%
overrepresentation of the mutated motif (as compared to what is likely by chance). For this reason, enrichment
is considered significant if it is greater than one. Motifs with enrichments <= 1 are not reported by MutaGene.

=============================
Examples
=============================

Option 1 - Search for motifs in mutational data

1. To search for all pre-identified motifs in sample1.maf using hg19

    $ mutagene motif search -i sample1.maf -g hg19 -s "="

    Sample Output and Interpretation:

    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    | sample                       | name       | motif      | strand | enrichment        | pvalue                 | mutations_low_est | mutations_high_est |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    | TCGA-50-6593-01A-11D-1753-08 | C>T in CpG | [C>T]G     | '='    | 4.586718025481874 | 1.0181609110804669e-06 | 15                | 18.0               |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+

    File "sample1.maf" contains one sample "TCGA-50-6593-01A-11D-1753-08"; from this sample 15-18
    mutations are estimated to be significantly contributed by the mutagenic process(es) involving C>T mutations in CpG motif ([C>T]G).
    The measures of significance are the enrichment and Fisher's Exact test pvalue calculations, where 0.05 is the threshold for statistical significance.

2. To search for the presence of the C[A>T] motif in sample1.maf using hg19

    $ mutagene motif search -i sample1.maf -g hg19 -m 'C[A>T]'

    No significant motif matches are found in the data, so nothing is reported.

3. To search sample1.maf for all preidentified motifs in mutagene on the transcription using hg19 and a window size of 20

    $ mutagene motif search -i sample1.maf -g hg19 -w 20 -s "+"

    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    | sample                       | name       | motif      | strand | enrichment        | pvalue                 | mutations_low_est | mutations_high_est |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    | TCGA-50-6593-01A-11D-1753-08 | APOBEC3G   | C[C>K]R    | '+'    |2.0770855332629354 | 0.022262032545564452   | 8                 | 14.0               |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    |TCGA-50-6593-01A-11D-1753-08  | C>T in CpG | [C>T]G     | '+'    |2.8697340043134436 |0.008360472489313148    | 7                 | 10.0               |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+

    File "sample1.maf" contains one sample "TCGA-50-6593-01A-11D-1753-08"; from this sample 8-14 mutations are estimated to be significantly contributed by the mutagenic process(es)
    involving APOBEC3G, where K represents the DNA bases G/T, and R represents the DNA bases A/G. 7-10 mutations are estimated to be significantly contributed by the mutagenic process(es) involving C>T mutations in CpG motif ([C>T]G).
    The measures of significance are the enrichment and Fisher's Exact test pvalue calculations, where 0.05 is the threshold for statistical significance.

To search sample2.vcf for all preidentified motifs in mutagene using hg19, searching for each of the motifs on the transcribed strand, non-trasncribed strand, plus both strands, and using a window size of plus/minus 30
bases from each mutation

    $ mutagene motif search -i sample2.vcf -g hg19 -w 30 -s "+-="

    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    | sample                       | name       | motif      | strand | enrichment        | pvalue                 | mutations_low_est | mutations_high_est |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    | VCF                          | APOBEC3G   | C[C>K]R    | '+'    |1.5208626215334309 | 9.767297094310342e-33  | 377               | 1099.0             |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    | VCF                          | APOBEC3G   | C[C>K]R    | '-'    |1.6115330339196352 |3.0535714666534214e-44  | 453               | 1193.0             |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    | VCF                          | APOBEC3G   | C[C>K]R    | '='    | 1.5665360537218949| 1.1734904382884064e-74 | 829               | 2292.0             |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    | VCF                          | C>T in CpG | [C>T]G     | '+'    |7.274092147503702  |0.0                     | 2029              | 2352.0             |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    | VCF                          | C>T in CpG | [C>T]G     | '-'    |4.248138083459255  |0.0                     | 1881              | 2460.0             |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    | VCF                          | C>T in CpG | [C>T]G     | '='    |11.074711617658798 |0.0                     | 4371              | 4804.0             |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    | VCF                          | Poly Eta   | W[A>T]     | '+'    |1.245342448790026  |0.013059702828698476    | 39                | 194.0              |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+
    | VCF                          | Poly Eta   | W[A>T]     | '='    |1.141805328027515  |0.020545858842258347    | 48                | 383.0              |
    +------------------------------+------------+------------+--------+-------------------+------------------------+-------------------+--------------------+

    File sample2.vcf was searched for all pre-identified motifs in MutaGene. Of these motifs, APOBEC3G and C>T in CpG
    were significantly present on the transcribed strand, non-transcribed strand, and both strands together.
    The presence of the Poly Eta motif was not significant on the non-transcribed strand but was significant on the
    transcribed stand and both the transcribed and non-transcribed strands together.

Option 2 - List all pre-identified motifs in MutaGene

    $ mutagene motif list

    The names and symbols for all pre-identified motifs in MutaGene will be listed
