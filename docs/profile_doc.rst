==================================
Description
==================================



==================================
Profile Subpackage Documentation
==================================

Command: mutagene profile calculate [--infile][--genome][--outfile][--help]

Required arguments (must be specified):

--infile, -i
    Input file in MAF or VCF format with one or multiple samples

--genome, -g
    Location of genome assembly file in 2bit format

Optional arguments (can be specified):

--outfile, -o
    Name of output file, will be generated  in TSV format

--help, -h
    show this help message and exit

===============================
How to Interpret Profile Output
===============================

The output will be in two columns.

The first column contains 96 context-dependent, mutation types (6 substitution types in 16 5'3' contexts).
The second column contains the total number of mutations that match the context-dependent,
mutation type in the left column.

Sample Output and Interpretation (from running mutagene profile calculate -i sample2.vcf -g hg19):

For the sake of conciseness, only a fraction of the total output is shown.

 - A[C>A]A 119
 - A[C>G]A 174
 - A[C>T]A 475
 - A[C>T]A 53

In sample2.vcf, 119 mutations match A[C>A]A, 174 match A[C>G]A, 475 match A[C>T]A, and 53 match A[C>T]A.