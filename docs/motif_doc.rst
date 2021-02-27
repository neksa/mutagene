===============================================================
Motif: Search for the presence of mutational motifs in samples
===============================================================

===============
1. Description
===============

------------------
1.1. Functionality
------------------

Use "mutagene motif" to search for the presence of mutational motifs in mutational data.

----------------------
1.2. Motif Definition
----------------------

A motif is defined as a characteristic pattern of a DNA mutation and its local DNA context. It is often associated with a specific carcinogen or a biological process.

--------------------------
1.3. Motif Representation
--------------------------

MutaGene represents motifs as strings of characters, where characters in brackets represent the single-base substitutions and characters outside brackets represent the unmutated DNA context. The motif must be in quotes to be recognized by MutaGene.

-------------------
1.4. Motif Examples
-------------------

"A[C>A]G" represents the DNA sequence "ACG" mutated into the DNA sequence "AAG"

"C[G>T]" represents the DNA sequence "CG" mutated into the DNA sequence "CT"

"[C>A]C" represents the DNA sequence "CC" mutated into the DNA sequence "AC"

------------------------------
1.5. Further Reading on Motifs
------------------------------

The publication `Mutational signatures and mutable motifs in cancer genomes <https://doi.org/10.1093/bib/bbx049>` describes motifs and their use.

*Note: if you installed MutaGene in a virtual environment, make sure you activate the virtual environment first.*

============================
2. Motif search command line
============================

To use the motif command, type 

``$ mutagene motif [arguments]``

You can always find help on the required arguments using the following command:

``$ mutagene motif -h``

============
3. Arguments
============

----------------------
3.1 Required arguments
----------------------

=========================   ============================================================  ====================
Argument                    Description                                                   Example
=========================   ============================================================  ====================
--infile INFILE             Input file in MAF or VCF format with 1 or multiple samples    --infile sample1.maf
                            (where INFILE is the sample filename with extension)
-i INFILE                   Short form of ``--infile INFILE`` argument                    -i sample1.maf 
--genome GENOME             Location of genome assembly file in 2bit format               --genome hg38.2bit   
                            (where GENOME is the filename)                    
-g GENOME                   Short form of ``--genome GENOME`` argument                    -g hg38.2bit 
=========================   ============================================================  ====================                                                                                                                                          

-----------------------------------------
3.2 Optional Arguments (can be specified)
-----------------------------------------

==========================  ============================================================  ============================
Argument                    Description                                                   Example
==========================  ============================================================  ============================
--motif MOTIF               Motif to search for, use the 'R[C>T]GY' syntax for the        --motif 'C[A>T]'
                            motif. Use quotes.
-m MOTIF                    Short form of ``--motif MOTIF``                               -m 'C[A>T]'
--input-format FORMAT       Input format, either MAF (the default) or VCF                 --input-format 'VCF'
-f FORMAT                   Short form of ``--input-format`` argument                     -i 'MAF'
--outfile OUTFILE           Name of output file, will be generated in TSV format          --outfile ../../out/out.tsv
                            (if this argument is not included output is to screen)
-o OUTFILE                  Short form of ``--outfile OUTFILE`` argument                  -o ../../out/out.tsv
--window-size WINDOW_SIZE   DNA local context window size for motif search                --window-size 30
                            (default setting is 50)\ :sup:`1`
-w WINDOW_SIZE              Short form of ``--window-size WINDOW_SIZE``                   -w 30
--strand {}                 Transcribed strand ("T"), non-transcribed ("N"), any ("A"),   --strand "T"
                            or all ("TNA" default)\ :sup:`2`
-s {}                       Short form of ``--strand {T,N,A,TNA}``                        -s "A"
--threshold THRESHOLD       Significance threshold for qvalues, default value=0.05        --threshold 0.03
-t THRESHOLD                Short form of ``--threshold THRESHOLD``                       -t 0.01
--save-motif-matches        Save mutations matching motif to a BED file                   --save-motif-matches out.bed
--test                      Statistical test to use, either Fisher (the default) or Chi2  --test 'Chi2'
==========================  ============================================================  ============================

1. Window Size Parameter Explanation: window size is defined as the number of DNA bases upstream and downstream from the mutated site not including the mutated site. Therefore, the effective length of the DNA sequence searched is 2 * window-size + 1. It may be advantageous to use a window size longer than the default, 50 bases, if the motif is longer than three nucleotides. Similarly, if the motif is shorter than three nucleotides, it may be advantageous to use a window size shorter than the default 50 bases. 

2. Strand Parameter Explanation: MutaGene can search for the presence of a motif on the transcribed or non-transcribed DNA strands or both strands. This information is gathered from the input file provided by the user. Analyzing the presence of a motif on a transcribed or non-transcribed strand is advantageous when a mutational process is known to have transcriptional strand bias. For instance, the APOBEC1/3A/B family is known to be associated with mutational processes that have a transcriptional strand bias of mutations. The transcription strand refers to the coding DNA strand, and the non-transcription strand refers to the template DNA strand.

=================================
4. Interpretation of Motif Output
=================================

If there no motifs significantly overrepresented in a predefined window, the output will say: "WARNING No significant motif matches found".

If the presence of a motif is significant in the data, the output will show a table with the following headers:

=============  ========================================================================================================
Header         Description
=============  ========================================================================================================
Sample         Name of Sample. If input file contains multiple samples, output will be stratified per sample.
Name           Name of motif. If -m/--motif argument is given, name will be "Custom motif".
Motif          Motif searched for in data
Strand         DNA Strand that motif was searched for on. "T": transcribed strand, "N": non-transcribed strand, "A":
               any strand, "TN": all strands.
Enrichment     Measure of motif's prevalence\ :sup:`3`
mut_low_est    Conservative estimate for a number of mutations (of total number in input file, mutational burden) that
               match a given motif
mut_high_est   Maximum number of mutations (of total number in input file) that match the motif
pvalue         Fisher's p-value for motif significance
qvalue         Fisher's p-value with Benjamini-Hochberg correction for multiple testing
=============  ========================================================================================================

3. How to Interpret Enrichment Output: Enrichment is modeled off of a risk ratio, meaning that a motif's enrichment is essentially a ratio between the probability of a motif appearing in a sample's DNA mutations and the probability of a motif appearing in a sample's DNA context. Enrichment minus one is equal to percent overrepresentation of a motif. For example, if enrichment is 1.5, it means that there is a 50% overrepresentation of the mutated motif. Motifs with enrichment values <= 1 are not reported by MutaGene.

===========
5. Examples
===========

-----------------------------------------------------------------------------------------
*5.1. Search for all pre-identfied motifs in sample1.maf using genome hg19 in any strand*
-----------------------------------------------------------------------------------------

---------------
5.1.1. Command
---------------

``$ mutagene motif -i sample1.maf -g hg19 -s "A"``

-------------------
5.1.2. Motif Output
-------------------

============================  ==========   ======  ==========  ==========  =======  =======  ==========  ==========  =========
sample                        mutagen      motif   strand      enrichment  mut_min  mut_max  odds_ratio  pvalue      qvalue   
============================  ==========   ======  ==========  ==========  =======  =======  ==========  ==========  =========
TCGA-50-6593-01A-11D-1753-08  C>T in CpG   [C>T]G  any strand  2.11727     10       18       2.42666     0.00169896  0.0118927
============================  ==========   ======  ==========  ==========  =======  =======  ==========  ==========  =========

--------------------------------
5.1.3. Interpretation of output
--------------------------------

File "sample1.maf" contains one sample "TCGA-50-6593-01A-11D-1753-08"; from this sample 10-18 mutations are estimated to be significantly contributed by the mutagenic process(es) involving C>T mutations in CpG motif ([C>T]G). The measures of significance are the enrichment and Fisher's Exact test pvalue calculations, where 0.05 is the threshold for statistical significance.

-----------------------------------------------------------------------------
*5.2. Search for the presence of the C[A>T] motif in sample1.maf using hg19*
-----------------------------------------------------------------------------

--------------
5.2.1. Command
--------------

``$ mutagene motif -i sample1.maf -g hg19 -m 'C[A>T]'``

-------------------
5.2.2. Motif Output
-------------------

No significant motif matches are found in the data, so nothing is reported.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
*5.3. Search sample2.vcf for all preidentified motifs in mutagene using hg19, searching for each of the motifs on the transcribed strand, non-transcribed strand, plus both strands, and using a window size of plus/minus 30 bases from each mutation*
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

--------------
5.3.1. Command
--------------

``$ mutagene motif -i sample2.vcf -g hg19 -w 30 -s "TNA" -f 'VCF'``

-------------------
5.3.2. Motif Output
-------------------

======  ==========  =======   ===============  ==========  =======  =======  ==========  ===========  ===========
sample  mutagen     motif     strand           enrichment  mut_min  mut_max  odds_ratio  pvalue       qvalue   
======  ==========  =======   ===============  ==========  =======  =======  ==========  ===========  ===========
VCF     C>T in CpG  [C>T]G    non-transcribed  4.24964     1882     2460     7.50456     0            0
VCF     C>T in CpG  [C>T]G    any strand       4.21114     3670     4812     7.38992     0            0
VCF     C>T in CpG  [C>T]G    transcribed      4.17202     1789     2352     7.27413     0            0
VCF     APOBEC3G    C[C>K]R   any strand       1.45802     721      2292     1.56628     1.17349e-74  6.16082e-74
VCF     APOBEC3G    C[C>K]R   non-transcribed  1.49104     393      1193     1.61104     3.05357e-44  1.2825e-43
VCF     APOBEC3G    C[C>K]R   transcribed      1.42365     328      1099     1.52034     9.7673e-33   3.41855e-32
VCF     Pol Eta     W[A>T]    transcribed      1.13727     24       194      1.24465     0.0130597    0.0391791
======  ==========  =======   ===============  ==========  =======  =======  ==========  ===========  ===========

-------------------------------
5.3.3. Interpretation of output
-------------------------------

File sample2.vcf was searched for all pre-identified motifs in MutaGene. Of these motifs, APOBEC3G and C>T in CpG were significantly present on the transcribed strand, non-transcribed strand, and both strands together. The presence of the Poly Eta motif was not significant on the non-transcribed strand but was significant on the transcribed stand and both the transcribed and non-transcribed strands together.
