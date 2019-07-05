==============================
What are motifs?
==============================

A motif is a characteristic pattern of DNA change, often associated with a carcinogen or a biological process. 
A motif's pattern includes a characteristic DNA context and a single-base substitution. 

----------------------------------------
How does MutaGene represent motifs?
----------------------------------------

MutaGene represents motifs in a string of characters, where the characters not in brackets represent the unmutated DNA context,
and the characters in brackets represent the single-base substitution.

Examples
---------- 

"A[C>A]G" represents the DNA sequence "ACG" mutated into the DNA sequence "AAG"

"C[G>T]" represents the DNA sequence "CG" mutated into the DNA sequence "CT"

"[C>A]C" represents the DNA sequence "CC" mutated into the DNA sequence "AC"

Motif Visualization
--------------------
.. image:: motif_viz.jpg
   :align: center
   :width: 100
   :height: 100
   :scale: 10 

--------------------------------------------------
Where can I read more about motifs?
--------------------------------------------------

This `publication <https://doi.org/10.1093/bib/bbx049>`_ describes motifs and their uses

==============================
Motif Subpackage Documentation
==============================

Option 1 - Search for motifs in genomic data

Usage: mutagene motif search [--infile][--genome][--motif][--outfile][--window-size][--strand][--help]

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

Option 2 - List all pre-identified motifs in MutaGene

Usage: mutagene motif list

--------------------------------------------------
Example Command-Line Arguments Using
--------------------------------------------------
