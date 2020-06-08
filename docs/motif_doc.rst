===============================================================
Motif: Search for the presence of mutational motifs in samples
===============================================================

===============
1. Description
===============

--------------------
1.1. Functionaility
--------------------

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
The publication `Mutational signatures and mutable motifs in cancer genomes <https://doi.org/10.1093/bib/bbx049>` describes motifs and their uses.

*Note: if you installed MutaGene in a virtual environment, make sure you activate the virtual environment first.*

-------------------
2. Motif command
-------------------

To use the motif command, type 

``$mutagene motif <action (search or list)>``

If search is specified, infile and genome are also required:

``$ mutagene motif search [arguments]``

You can always find help on the required arguments using the following command:

``$ mutagene motif search -h``
or
``$mutagene motif list -h``

------------
3. Arguments
------------

**3.1.Command:** ``$mutagene motif <action (search or list)>``

followed by the required arguments from the command line. 

**3.2.Required Arguments (must be specified):**

Motif function requires:
``$mutagene motif <action (search or list)>``
If search is specified, infile and genome are also required

=========================   ============================================================  ====================
Argument                    Description                                                   Example
=========================   ============================================================  ====================
--infile INFILE             Input file in MAF or VCF format with 1 or multiple samples     --infile sample1.maf
                            (where INFILE is the sample filename with extension)
-i INFILE                   Short form of --infile INFILE argument                         -i sample1.maf 
--genome GENOME             Location of genome assembly file in 2bit format                --genome hg38.2bit   
                            (where GENOME is the filename)                    
-g GENOME                   Short form of --genome GENOME argument                         -g hg38.2bit                      
=========================   ============================================================  ====================                                                                                                                                          


**3.3.Optional Arguments (can be specified):**

==========================  =============================================================  ============================
Argument                    Description                                                    Example
==========================  =============================================================  ============================
--outfile [OUTFILE]         Name of output file, will be generated in TSV format            --outfile ../../out/out.tsv
                            (if this argument is not included output is to screen)
-o [OUTFILE]                Short form of --outfile [OUTFILE] argument                      -o ../../out/out.tsv
--motif MOTIF               Motif to search for, use the 'R[C>T]GY syntax for the           --motif 'C[A>T]'
                            motif. Use quotes
-m MOTIF                    Short form of --motif MOTIF                                     -m 'C[A>T]'
window-size WINDOW_SIZE     Context window size for motif search                            --window-size 30
                            (default setting is 50)\ :sup:`1`
-w WINDOW_SIZE              Short form of window-size WINDOW_SIZE                           -w 30
--strand {+,-,=,+-=}        Transcribed strand (+), non-transcribed (-), any (=),           --strand "+-="
                            or all (+-= default)\ :sup:`2`
-s {+,-,=,+-=}              Short form of --strand {+,-,=,+-=}                              -s "+-="
==========================  =============================================================  ============================

1. Window Size Parameter Explanation: MutaGene counts window size as the number of DNA bases searched for from the first base of the DNA sequence gathered up to but not including the mutated base. Therefore, the effective length of the DNA sequence searched is 2 * window-size + 1. It may be advantageous to use a window size longer than the default 50 bases if the motif is longer than three nucleotides,
as this motif is likely to appear less frequently in the DNA context. Similarly, if the motif is shorter than three nucleotides,
it may be advantageous to use a window size shorter than the default 50 bases, as the motif is likely to appear in DNA more frequently.
2. Strand Parameter Explanation: MutaGene can search for the presence of a motif on the transcribed or non-transcribed DNA strands or both strands. This information is gathered from the input file provided by the user. Analyzing for the presence on a transcribed or non-transcribed strand is advantageous when a mutational process is known to have mutations with a transcriptional strand bias. For instance, the APOBEC1/3A/B family is known to be associated with mutational processes that have a transcriptional strand bias of mutations in exons. The transcription strand refers to the coding DNA strand, and the non-transcription strand refers to the template DNA strand.

=======

The publication `Mutational signatures and mutable motifs in cancer genomes <https://doi.org/10.1093/bib/bbx049>` describes motifs and their use.

*Note: if you installed MutaGene in a virtual environment, make sure you activate the virtual environment first.*

-------------------
2. Motif search command line
-------------------

To use the motif command, type 

``$ mutagene motif [arguments]``

You can always find help on the required arguments using the following command:

``$ mutagene motif -h``

------------
3. Arguments
------------

3.1 Required arguments
=========================   ============================================================  ====================
Argument                    Description                                                   Example
=========================   ============================================================  ====================
--infile INFILE             Input file in MAF or VCF format with 1 or multiple samples     --infile sample1.maf
                            (where INFILE is the sample filename with extension)
-i INFILE                   Short form of --infile INFILE argument                         -i sample1.maf 
--genome GENOME             Location of genome assembly file in 2bit format                --genome hg38.2bit   
                            (where GENOME is the filename)                    
-g GENOME                   Short form of --genome GENOME argument                         -g hg38.2bit 
--motif MOTIF               Motif to search for, use the 'R[C>T]GY syntax for the           --motif 'C[A>T]'
                            motif. Use quotes
-m MOTIF                    Short form of --motif MOTIF                                     -m 'C[A>T]'
=========================   ============================================================  ====================                                                                                                                                          


3.2 Optional Arguments (can be specified)

==========================  =============================================================  ============================
Argument                    Description                                                    Example
==========================  =============================================================  ============================
--outfile [OUTFILE]         Name of output file, will be generated in TSV format            --outfile ../../out/out.tsv
                            (if this argument is not included output is to screen)
-o [OUTFILE]                Short form of --outfile [OUTFILE] argument                      -o ../../out/out.tsv
window-size WINDOW_SIZE     DNA local context window size for motif search                  --window-size 30
                            (default setting is 50)\ :sup:`1`
-w WINDOW_SIZE              Short form of window-size WINDOW_SIZE                           -w 30
--strand {}                 Transcribed strand ("T"), non-transcribed ("N"), any ("A"),     --strand ""
                            or all ("NT" default)\ :sup:`2`
-s {}                       Short form of --strand {+,-,=,+-=}                              -s ""
--save-motif-matches        Save mutations matching motif to a BED file
SAVE_MOTIF_MATCHES
                        
==========================  =============================================================  ============================

1. Window Size Parameter Explanation: window size is defined as the number of DNA bases upstream and downstream from the mutated site not including the mutated site. Therefore, the effective length of the DNA sequence searched is 2 * window-size + 1. It may be advantageous to use a window size longer than the default, 50 bases, if the motif is longer than three nucleotides. Similarly, if the motif is shorter than three nucleotides, it may be advantageous to use a window size shorter than the default 50 bases. 

2. Strand Parameter Explanation: MutaGene can search for the presence of a motif on the transcribed or non-transcribed DNA strands or both strands. This information is gathered from the input file provided by the user. Analyzing the presence of a motif on a transcribed or non-transcribed strand is advantageous when a mutational process is known to have transcriptional strand bias. For instance, the APOBEC1/3A/B family is known to be associated with mutational processes that have a transcriptional strand bias of mutations. The transcription strand refers to the coding DNA strand, and the non-transcription strand refers to the template DNA strand.

---------------------------------
4. Interpretation of Motif Output
---------------------------------
If there no motifs significantly overrepresented in a predefined window, the output will say: "WARNING No significant motif matches found".

If the presence of a motif is significant in the data, the output will show a table with the following headers:

=======================================================================================================================
Header         Description
=======================================================================================================================
Sample         Name of Sample. If input file contains multiple samples, output will be stratified per sample.
Name           Name of motif. If -m/--motif argument is given, name will be "Custom motif".
Motif          Motif searched for in data
Strand         DNA Strand that motif was searched for on. 'T': transcribed strand, 'N': non-transcribed strand, "A": any                      strand, "TN": all strands.
Enrichment     Measure of motif's prevalence \ :sup:`a`
mut_low_est    Conservative estimate for a number of mutations (of total number in input file, mutational burden) that match a                given motif
mut_high_est   Maximum number of mutations (of total number in input file) that match the motif
pvalue         Fisher's p-value for motif significance
qvalue         Fisher's p-value with Benjamini-Hochberg correction for multiple testing
=============  =======================================================================================================================

a. How to Interpret Enrichment Output: Enrichment is modeled off of a risk ratio, meaning that a motif’s enrichment is essentially a ratio between the probability of a motif appearing in a sample’s DNA mutations and the probability of a motif appearing in a sample’s DNA context. Enrichment minus one is equal to percent overrepresentation of a motif. For example, if enrichment is 1.5, it means that there is a 50% overrepresentation of the mutated motif. Motifs with enrichment values <= 1 are not reported by MutaGene.

-----------
5. Examples
-----------------------------------------------------------------------------------------
*5.1. Search for all pre-identfied motifs in sample1.maf using genome hg19 in any strand*
-----------------------------------------------------------------------------------------

---------------
5.1.1. Command
---------------

``$ mutagene motif -i sample1.maf -g hg19 -s "="``

-------------------
5.1.2. Motif Output
-------------------

============================  ===========  ======  ======  =================  ======================  ===========  ============
sample                        name         motif   strand  enrichment         pvalue                  mut_low_est  mut_high_est   
============================  ===========  ======  ======  =================  ======================  ===========  ============
TCGA-50-6593-01A-11D-1753-08  C>T in CpG   [C>T]G  '='     4.586718025481874  1.0181609110804669e-06  15           18.0
============================  ===========  ======  ======  =================  ======================  ===========  ============ 

--------------------------------
5.1.3. Interpretation of output
--------------------------------
File "sample1.maf" contains one sample "TCGA-50-6593-01A-11D-1753-08"; from this sample 15-18 mutations are estimated to be significantly contributed by the mutagenic process(es) involving C>T mutations in CpG motif ([C>T]G). The measures of significance are the enrichment and Fisher's Exact test pvalue calculations, where 0.05 is the threshold for statistical significance.

-----------------------------------------------------------------------------
*5.2. Search for the presence of the C[A>T] motif in sample1.maf using hg19*
-----------------------------------------------------------------------------

-------------
5.2.1. Command
-------------

``$ mutagene motif -i sample1.maf -g hg19 -m 'C[A>T]'``

-------------------
5.2.2. Motif Output
-------------------

No significant motif matches are found in the data, so nothing is reported.

--------------------------------------------------------------------------------------------------------------------------
*5.3. Search sample1.maf for all preidentified motifs in mutagene on the transcription using hg19 and a window size of 20*
--------------------------------------------------------------------------------------------------------------------------

-------------
5.3.1. Command
-------------

``$ mutagene motif -i sample1.maf -g hg19 -w 20 -s "+"``

-------------------
5.3.2. Motif Output
-------------------

============================  ===========  =======  ======  =================  ======================  ===========  ============
sample                        name         motif    strand  enrichment         pvalue                  mut_low_est  mut_high_est   
============================  ===========  =======  ======  =================  ======================  ===========  ============
TCGA-50-6593-01A-11D-1753-08  APOBEC3G     C[C>K]R  '+'    2.0770855332629354  0.022262032545564452   8             14.0
TCGA-50-6593-01A-11D-1753-08  C>T in CpG   [C>T]G   '+'    2.8697340043134436  0.008360472489313148   7             10.0
============================  ===========  =======  ======  =================  ======================  ===========  ============

--------------------------------
5.3.3. Interpretation of output
--------------------------------

File "sample1.maf" contains one sample "TCGA-50-6593-01A-11D-1753-08"; from this sample 8-14 mutations are estimated to be significantly contributed by the mutagenic process(es) involving APOBEC3G, where K represents the DNA bases G/T, and R represents the DNA bases A/G. 7-10 mutations are estimated to be significantly contributed by the mutagenic process(es) involving C>T mutations in CpG motif ([C>T]G).
The measures of significance are the enrichment and Fisher's Exact test pvalue calculations, where 0.05 is the threshold for statistical significance.


-----------------------------------------------------------------------------
*5.4. Search sample2.vcf for all preidentified motifs in mutagene using hg19, searching for each of the motifs on the transcribed strand, non-transcribed strand, plus both strands, and using a window size of plus/minus 30
bases from each mutation*
-----------------------------------------------------------------------------

-------------
5.4.1. Command
-------------

``$ mutagene motif -i sample2.vcf -g hg19 -w 30 -s "+-="``

-------------------
5.4.2. Motif Output
-------------------

======  ===========  =======  ======  ==================  ======================  ===========  ============
sample  name         motif    strand  enrichment          pvalue                  mut_low_est  mut_high_est   
======  ===========  =======  ======  ==================  ======================  ===========  ============
VCF     APOBEC3G     C[C>K]R  '+'     1.5208626215334309  9.767297094310342e-33   377          1099.0
VCF     APOBEC3G     C[C>K]R  '-'     1.6115330339196352  3.0535714666534214e-44  453          1193.0
VCF     APOBEC3G     C[C>K]R  '='     1.5665360537218949  1.1734904382884064e-74  829          2292.0
VCF     C>T in CpG   [C>T]G   '+'     7.274092147503702   0.0                     2029         2352.0
VCF     C>T in CpG   [C>T]G   '-'     4.248138083459255   0.0                     1881         2460.0
VCF     C>T in CpG   [C>T]G   '='     11.074711617658798  0.0                     4371         4804.0
VCF     Poly Eta     W[A>T]   '+'     1.245342448790026   0.013059702828698476    39           194.0
VCF     Poly Eta     W[A>T]   '='     1.141805328027515   0.020545858842258347    48           383.0
======  ===========  =======  ======  ==================  ======================  ===========  ============

--------------------------------
5.4.3. Interpretation of output
--------------------------------

File sample2.vcf was searched for all pre-identified motifs in MutaGene. Of these motifs, APOBEC3G and C>T in CpG were significantly present on the transcribed strand, non-transcribed strand, and both strands together. The presence of the Poly Eta motif was not significant on the non-transcribed strand but was significant on the transcribed stand and both the transcribed and non-transcribed strands together.

-----------------------------------------------------------------------------------------
*5.5. List all pre-identified motifs in MutaGene*
-----------------------------------------------------------------------------------------

-------------
5.5.1. Command
-------------

``$ mutagene motif list``

-------------------
5.5.2. Motif Output
-------------------

============  =======
name          symbol
============  =======
APOBEC1/3A/B  T[C>K]W
APOBEC3G      C[C>K]R
C>T in CpG    [C>T]G
UV Light      Y[C>T]
Pol Eta       W[A>T]
AID           W[R>S]C
============  =======

--------------------------------
5.5.3. Interpretation of output
--------------------------------

The names and symbols for all pre-identified motifs in MutaGene are listed.


=======

--------------------------------
5.4.3. Interpretation of output
--------------------------------

File sample2.vcf was searched for all pre-identified motifs in MutaGene. Of these motifs, APOBEC3G and C>T in CpG were significantly present on the transcribed strand, non-transcribed strand, and both strands together. The presence of the Poly Eta motif was not significant on the non-transcribed strand but was significant on the transcribed stand and both the transcribed and non-transcribed strands together.
