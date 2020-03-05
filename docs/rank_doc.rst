=====================================================
Rank: Identifying potential drivers in samples
=====================================================
-----------
Description
-----------
Use mutagene rank to analyze genes and compare observed mutational frequencies to expected background mutability to identify potential drivers. To read more about the MutaGene's ranking method, read the publication
`Finding driver mutations in cancer: Elucidating the role of background mutational processes <https://doi.org/10.1371/journal.pcbi.1006981>`_

*Note: if you installed MutaGene in a virtual environment, make sure you activate the virtual environment first.*

-------------------
1. Rank command
-------------------

To use the rank command, type 

``$ mutagene rank``

followed by the required arguments from the command line. You can always find help on the required arguments using the following command:

``$ mutagene rank -h``

------------
2. Arguments
------------

**2.1.Command:** mutagene rank [arguments]

**2.2.Required Arguments (must be specified):**

=========================   ============================================================  ====================
Argument                    Description                                                   Example
=========================   ============================================================  ====================
--infile INFILE             Input file in MAF format                                       --infile sample1.maf
                            (where INFILE is the sample filename with extension)
-i INFILE                   Short form of --infile INFILE argument                         -i sample1.maf 
--genome GENOME             Location of genome assembly file in 2bit format                --genome hg19
-g GENOME                   Short form of --genome GENOME argument                         -g hg19
--cohort COHORT             Name of cohort with observed mutations                         --cohort gcb_lymphomas
-c COHORT                   Short form of --cohort COHORT argument                         -c gcb_lymphomas
=========================   ============================================================  ====================                                                                                                                                   

**2.3.Optional Arguments (can be specified):**

=========================================  =================================================  ==================================
Argument                                   Description                                        Example
=========================================  =================================================  ==================================
--outfile OUTFILE                          Name of output file, will be generated in           --outfile out.tsv
                                           TSV format  (If this argument is not included,
                                           output is to screen)   
-o OUTFILE                                 Short form of --outfile OUTFILE                     -o out.tsv
--profile PROFILE                          Override profile to calculate mutability
                                           (may also describe cohort size)
-p PROFILE                                 Short form of --profile PROFILE
--nsamples NSAMPLES                        Override cohort size    
-n NSAMPLES                                Short form of --nsamples
--threshold-driver THRESHOLD_DRIVER        BScore threshold between Driver and Pontential 
                                           Driver mutations
-td THRESHOLD_DRIVER                       Short form of --threshold-driver
--threshold-passenger THRESHOLD_PASSENGER  BScore threshold between Potential Driver and 
                                           Passenger mutations
-tp THRESHOLD_PASSENGER                    Short form of --threshold-passenger
--cohorts-file COHORTS_FILE                Location of tar.gz container or directory 
                                           for cohorts
=========================================  =================================================  ==================================  

--------------------------------
3. Interpretation of Rank Output
--------------------------------

The output will show a table with the headers described in table below. 

===================  =======================================================================================================
Output Table Header  Description    
===================  =======================================================================================================
gene                 Name of gene with mutation
mutation             Expressed as, eg. Y99F, ie. amino acid tyrosine (Y) replaced by phenylalanine (F) at position 99  
mutability           Expected mutation rate in a particular DNA context
observed             Observed mutational frequencies
bscore               A binomial p-value for the observed number of occurences of mutation in comparison to the expected
                     mutability that is defined by the local DNA context of the mutated nucleotide
qvalue               Bscore corrected for multiple testing with Benjamini-Hochberg FDR method
label                Prediction of cancer drivers, Potential drivers, and Passengers is based on the thresholds established
                     for the Bscore optimized using this benchmark datasets.
===================  =======================================================================================================

-----------
4. Examples
-----------

---------------------------------------------------------------------------------------------------
*4.1. Use mutagene rank to analyze genes in sample1.maf using genome hg19 and cohort gcb_lymphomas*
---------------------------------------------------------------------------------------------------

-------------
4.1.1.Command
-------------

``$ mutagene rank -i sample1.maf -g hg19 -c gcb_lymphomas``

-------------------------------------------------------
4.1.2.Rank Output (only first 5 results are shown here)
-------------------------------------------------------

========  =========  =======================  ========  =======================  =====================  ======    
gene      mutation   mutability               observed  bscore                   qvalue                 label   
========  =========  =======================  ========  =======================  =====================  ======  
BOD1L     T2810S     8.09229314668869e-08     1         3.6415254329818015e-06   5.577139539596271e-05  Driver
TEX15     V2686E     8.540363127806927e-08    1         3.843156186679522e-06    5.577139539596271e-05  Driver
GRINA     Y99F       8.540363127806927e-08    1         3.843156186679522e-06    5.577139539596271e-05  Driver
N4BP2L2   K143I      1.0351675938657934e-07   1         4.658243563849532e-06    5.577139539596271e-05  Driver
ZC3H3     R59G       1.1254702103613567e-07   1         5.06460340648271e-06     5.577139539596271e-05  Driver
========  =========  =======================  ========  =======================  =====================  ======   

------------------------------------------------------------------------------------------------------------------------------------
*4.2. Use mutagene rank to analyze genes in sample1.maf using genome hg19 and cohort gcb_lymphomas with a BScore threshold of 0.0003 between Potential Driver and Passenger mutations *
------------------------------------------------------------------------------------------------------------------------------------

-------------
4.2.1.Command
-------------

``$ mutagene rank -i sample1.maf -g hg19 -c gcb_lymphomas -tp 0.0003``

-----------------------------------------------------------------------------------
4.2.2.Rank Output (only 4 results around potential driver and passenger are shown here)
-----------------------------------------------------------------------------------

========  =========  =======================  ========  =======================  ======================  ================    
gene      mutation   mutability               observed  bscore                   qvalue                  label   
========  =========  =======================  ========  =======================  ======================  ================  
WNT8B     R231C      6.280123772905988e-06    1         0.00028256652774017057   0.00029419110008391177  Potential Driver
ATXN1     P109L      6.280123772905988e-06    1         0.00028256652774017057   0.00029419110008391177  Potential Driver
OR2T12    P180P      6.797840069627803e-06    1         0.0003058570590671567    0.0003096214536402909   Passenger
GPR77     S333S      6.797840069627803e-06    1         0.0003058570590671567    0.0003096214536402909   Passenger
========  =========  =======================  ========  =======================  ======================  ================

--------------------------------------------------------------------------------------------------------------------------------------
*4.3. Use mutagene rank to analyze genes in sample1.maf using genome hg19 and cohort gcb_lymphomas with a BScore threshold of 0.000009 between Driver and Potential Driver mutations *
--------------------------------------------------------------------------------------------------------------------------------------

-------------
4.2.1.Command
-------------

``$ mutagene rank -i sample1.maf -g hg19 -c gcb_lymphomas -td 0.000009``

-----------------------------------------------------------------------------------
4.2.2.Rank Output (only 4 results around driver and potential driver are shown here)
-----------------------------------------------------------------------------------

========  =========  =======================  ========  =======================  ======================  ================    
gene      mutation   mutability               observed  bscore                   qvalue                  label   
========  =========  =======================  ========  =======================  ======================  ================  
C1orf69   E244V      1.9422490304954465e-07   1         8.740083291253642e-06    5.577139539596271e-05   Driver
PARD3B    E1055V     1.9422490304954465e-07   1         8.740083291253642e-06    5.577139539596271e-05   Driver
KIF21B    L517V      2.1106070979826086e-07   1         9.497687839898163e-06    5.577139539596271e-05   Potential Driver
KIAA1409  R294L       2.1106070979826086e-07  1         9.497687839898163e-06    5.577139539596271e-05   Potential Driver
========  =========  =======================  ========  =======================  ======================  ================

----------------------------------------------------------------------------------------------------------------------------
*4.3. Use mutagene rank to analyze genes in sample1.maf using genome hg19 and cohort gcb_lymphomas with a cohort size of 20*
----------------------------------------------------------------------------------------------------------------------------

-------------
4.3.1.Command
-------------

``$ mutagene rank -i sample1.maf -g hg19 -c gcb_lymphomas -n 20``

-------------------------------------------------------
4.3.2.Rank Output (only first 5 results are shown here)
-------------------------------------------------------

========  =========  =======================  ========  =======================  =====================  ======    
gene      mutation   mutability               observed  bscore                   qvalue                 label   
========  =========  =======================  ========  =======================  =====================  ======  
BOD1L     T2810S     1.7803044916053778e-07   1         3.7386327764622237e-06   5.725863260405688e-05  Driver
TEX15     V2686E     1.8788798872293455e-07   1         3.945640349792222e-06    5.725863260405688e-05  Driver
GRINA     Y99F       1.8788798872293455e-07   1         3.945640349792222e-06    5.725863260405688e-05  Driver
N4BP2L2   K143I      2.2773687058386116e-07   1         4.782463390819526e-06    5.725863260405688e-05  Driver
ZC3H3     R59G       2.4760344619068064e-07   1         5.199659495456503e-06    5.725863260405688e-05  Driver
========  =========  =======================  ========  =======================  =====================  ======   
