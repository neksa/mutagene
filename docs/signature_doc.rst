==================================
Signature Subpackage Documentation
==================================

Command: mutagene signature identify [--signatures {5,10,30}][--infile][--genome][--outfile][--help]

Required arguments (must be specified):

--signatures, -s
    Collection of signatures to use, choose from 5, 10, or 30

--infile, -i
    Input file in MAF or VCF format with one or multiple samples

--genome, -g
    Location of genome assembly file in 2bit format

Optional arguments (can be specified):

--outfile, -o
    Name of output file, will be generated  in TSV format

--help, -h
    show this help message and exit

Examples:
    1. mutagene signature identify -s 5 -i sample1.maf -g hg19 --> identifies presence of MutaGene-5 signatures in
       sample1.maf using hg19

    2. mutagene signature identify -s 10 -i sample2.vcf -g hg18 ---> identifies presence of MutaGene-10 signatures in
       sample2.vcf using hg18

----------------------------------------------------
What Signatures Can I Analyze Using MutaGene Package
----------------------------------------------------

The MutaGene signature package allows for the analysis of 3 different "bundles" of mutational signatures: MutaGene-5 Signatures, MutaGene-10 Signatures, and Cosmic-30 Signatures.
MutaGene-5 contains 5 signatures, MutaGene-10 contains 10 signatures, and Cosmic-30 contains 30 signatures.

Here is the `link <https://www.ncbi.nlm.nih.gov/research/mutagene/signatures#mutational_signatures>`_ to see more about the 3 available signature bundles.

=================================
How to Interpret Signature Output
=================================

The output will be in two columns.

The first column will have whole numbers from 1-5 if you chose MutaGene-5, 1-10 if you chose
MutaGene-10, or 1-30 if you chose COSMIC-30. Each number is the signature number.

The second column will have decimal numbers; these decimals represent the contribution (aka the proportion of
mutations estimated to be contributed by the mutagenic process(es) represented by the signature.

Sample Output and Interpretation (from running mutagene signature identify -s 5 -i sample1.maf -g hg19):

+---+---------+
| 1 | 0.5361  |
+---+---------+
| 2 | 0.0556  |
+---+---------+
| 3 | 0.1105  |
+---+---------+
| 4 | 0.11183 |
+---+---------+
| 5 | 0.0506  |
+---+---------+

This output means that 53.61% of input mutations are estimated to be contributed by MutaGene-5 Signature 1,
5.56% of input mutations are estimated to be contributed by MutaGene-5 Signature 2, etc.