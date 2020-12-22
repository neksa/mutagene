=====================================================
Rank: identifying potential driver mutations
=====================================================
-----------
Description
-----------
"mutagene rank" module ranks mutations with respect to their driver statuses. The method requires three input parameters: the background mutability model for each nucleotide or codon, a number of samples where a given mutation was observed (mutational frequency) and the overall number of samples in a given cohort of patients. The background mutability model, mutational frequency and number of samples are specified and taken from the input file by default. Arguments below can overwrite the default.

Please cite the MutaGene ranking method as 
Anna-Leigh Brown, Minghui Li, Alexander Goncearenco, Anna R Panchenko
"Finding Driver Mutations in Cancer: Elucidating the Role of Background Mutational Processes" Plos Comp Biol 15 (4), e1006981
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006981

*Note: if you installed MutaGene in a virtual environment, make sure you activate the virtual environment first.*


-------------------
1. Rank 
-------------------

To rank mutations, type 

$ mutagene rank

followed by the required arguments, you can always find help on the required arguments using the following command:

$ mutagene rank -h

------------
2. Arguments
------------

**2.1.Command:** mutagene rank [arguments]

**2.2.Required Arguments (must be specified):**

=========================   ============================================================  ====================
Argument                    Description                                                   Example
=========================   ============================================================  ====================
--infile INFILE             Input file of mutations to be ranked in MAF format            --infile sample1.maf
                            (where INFILE is the sample(s) filename with extension)
-i INFILE                   Short form of --infile INFILE argument                         -i sample1.maf 
--genome GENOME             Location of genome assembly file in 2bit format                --genome hg19
-g GENOME                   Short form of --genome GENOME argument                         -g hg19

=========================   ============================================================  ====================                                                                                                                                   

**2.3.Optional Arguments (can be specified):**

=========================================  =================================================  ==================================
Argument                                   Description                                        Example
=========================================  =================================================  ==================================
--outfile OUTFILE                          Name of output file, will be generated in           --outfile out.tsv
                                           TSV format  (If this argument is not included,
                                           output is to screen)   
-o OUTFILE                                 Short form of --outfile OUTFILE                     -o out.tsv
--cohort COHORT                            Name of precalculated cohort which overwrites  
                                           input sample(s). If cohort is specified, all three  --cohort gcb_lymphomas
                                           method's input parameters will be derived from it.                                                                       
-c COHORT                                  Short form of --cohort COHORT argument              -c gcb_lymphomas
--profile PROFILE                          Overrides background mutability model                                          
-p PROFILE                                 Short form of --profile PROFILE
--nsamples NSAMPLES                        Overrides the number of samples in cohort          --nsamples 20
-n NSAMPLES                                Short form of --nsamples                           -n 20
--threshold-driver THRESHOLD_DRIVER        BScore threshold between Driver and Pontential     --threshold-driver 0.000009
                                           Driver mutations
-td THRESHOLD_DRIVER                       Short form of --threshold-driver                   -td 0.000009
--threshold-passenger THRESHOLD_PASSENGER  BScore threshold between Potential Driver and      --threshold-passenger 0.0003
                                           Passenger mutations
-tp THRESHOLD_PASSENGER                    Short form of --threshold-passenger                -tp 0.0003
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
mutability           Expected mutation rate in a given DNA context
observed             Number of cancer samples where this mutation was observed
bscore               A binomial p-value for the observed number of occurences of mutation in comparison to the expected
                     mutability that is defined by the local DNA context of the mutated nucleotide
qvalue               Bscore corrected for multiple testing with Benjamini-Hochberg FDR method
label                Prediction of cancer drivers, Potential drivers, and Passengers is based on the thresholds established
                     for the Bscore optimized using this benchmark datasets. This is a rather arbitrary threshold.
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
