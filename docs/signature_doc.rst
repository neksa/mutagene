=====================================================
Identify: Identifying mutational profiles in samples
=====================================================
-----------
Description
-----------
Use mutagene signature to search for the presence of mutational signatures in mutational data.

*Note: if you installed MutaGene in a virtual environment, make sure you activate the virtual environment first.*

-------------------
1. Identify command
-------------------

To use the identify command, type 

``$ mutagene signature identify``

followed by the required arguments from the command line. You can always find help on the required arguments using the following command:

``$ mutagene signature identify -h``

------------
2. Arguments
------------

**2.1.Command:** mutagene signature identify [arguments]

**2.2.Required Arguments (must be specified):**

=========================   ============================================================  ====================
Argument                    Description                                                   Example
=========================   ============================================================  ====================
--infile INFILE             Input file (with one or multiple samples) in VCF or MAF format  --infile PD3851a.vcf
                            (where INFILE is a sample filename)
-i INFILE                   Short form of --infile INFILE argument                         -i PD3851a.vcf 
--genome GENOME             Location of a genome assembly file in 2bit format                --genome hg38.2bit   
                            (where GENOME is the filename)                    
-g GENOME                   Short form of --genome GENOME argument                         -g hg38.2bit                      
--signatures {5,10,30,49}   Collection of mutational signatures (include MutaGene5 (5),       --signatures10 
                            MutaGene10(10), Cosmic30(30) and Cosmic49(49))\ :sup:`1`
-s {5,10,30,49}             Short form of --signatures {5,10,30,49}                        -s10  
=========================   ============================================================  ====================                                                                                                                                          

1. The MutaGene signature package allows for the analysis of mutational processes operating in a given sample or a cohort using four different sets of mutational signatures: MutaGene5, MutaGene10, Cosmic30 and Cosmic49 signature sets.
Read more about the MutaGene signature packag at https://www.ncbi.nlm.nih.gov/research/mutagene/signatures#mutational_signatures.
Read more about Cosmic signatures at https://cancer.sanger.ac.uk/cosmic/signatures/.


**2.3.Optional Arguments (can be specified):**

==========================  =============================================================  ============================
Argument                    Description                                                    Example
==========================  =============================================================  ============================
--input-format {MAF,VCF}    Input file format: MAF, VCF                                     --input-format VCF
                            (MAF is a default format)
-f {MAF,VCF}                Short form of --input-format {MAF,VCF} argument                 -f VCF 
--outfile [OUTFILE]         Name of output file, will be generated in TSV format            --outfile ../../out/out.tsv
                            
-o [OUTFILE]                Short form of --outfile [OUTFILE] argument                      -o ../../out/out.tsv
==========================  =============================================================  ============================


**2.4.Advanced arguments:**

==========================  ===================================================================  =========================
Argument                    Description                                                          Example
==========================  ===================================================================  =========================
--method [METHOD]           Method defines a target function minimized in the optimization 
                            procedure                                                            -method frobenius
                            (default method is mlez)\ :sup:`1`
-m [METHOD]                 Short form of --method [METHOD] argument                             -m frobenius
--no-unexplained-variance   not to account for unexplained variance from context-independent     --no-unexplained-variance
                            mutational processes and unknown signatures                                                       
-U                          Short form of --no-unexplained-variance argument                     -U
--bootstrap                 Use the bootstrap to calculate 95% confidence intervals              --bootstrap
                            for VCF files                                                        
-b                          Short form of --bootstrap argument                                   -b
==========================  ===================================================================  =========================

The bootstraping is done by random sampling with replacement of mutations from each sample repeated 100 times. The bootstrap estimates of the mean and 95% confidence intervals are reported.

1. Available methods

================= ===================================================================================== 
Argument           Method      
================= ===================================================================================== 
frobenius          Frobenius  
frobeniuszero      FrobeniusZero 
js                 js
divergencejs       divergencejs
mle                NegLogLik (Maxiumum Likelihood Estimation maximizes LogLik or minimizes NegLogLik) 
mlez               NegLogLik (MLE with added context-independent signatures)
compat             MegLogLikOld (MLE with added context-independent signatures (old compatability mode))
aicc               AICc (AIC corrected for small samples) 
bic                BIC
aiccz              BIC (BIC with added context-independent signatures)  
================= =====================================================================================


-----------
3. Examples
-----------
*3.1. Search for the presence of MutaGene-10 signatures in PD3851a.vcf using hg38 and default method mlez (method mle gives same output for this input)*
-------
3.1.1.Command
-------

``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10``

---------------
3.1.2.Identify Output
---------------

=======  ============  ============  =========== 
sample     signature     exposure    mutations   
=======  ============  ============  =========== 
VCF        2             0.0935      80
VCF        3             0.0392      33
VCF        4             0.0074      6
VCF        5             0.0728      62
VCF        6             0.1362      116
VCF        7             0.0118      10
VCF        8             0.0552      47
VCF        9             0.0271      23
VCF        10            0.0121      10
=======  ============  ============  =========== 


*3.2. Calculate the mutational profile for sample1.maf using -g hg38.2bit and MutaGene-5 signature set:*
-------
3.2.1.Command
-------

``$ mutagene signature identify -i sample1.maf -g hg38.2bit -s5``

---------------
3.2.2.Identify Output
---------------

=============================  ===========  ==========  ===========  
sample                          signature    exposure    mutations
=============================  ===========  ==========  ===========
TCGA-50-6593-01A-11D-1753-08    2            0.0348      5
TCGA-50-6593-01A-11D-1753-08    3            0.0691      11
=============================  ===========  ==========  ===========  

*3.3. Calculate the mutational profile for PD3851a.vcf using -g hg38.2bit (both in samples folder) and MutaGene-10 signature set and send output to a file out.tsv in out folder:*
-------
3.3.1.Command
-------

``$ mutagene signature identify -i ../../samples/PD3851a.vcf -g ../../samples/hg38.2bit -f VCF -s10 -o ../../out/out.tsv``

---------------
3.3.2.Identify Output
---------------
As for example (3.1) except sent to file out.tsv instead of screen.

*3.4. Calculate the mutational profile for PD3851a.vcf using -g hg38.2bit and MutaGene-10 signature set using the bootstrap to calculate confidence intervals:*
-------
3.4.1.Command
-------

``$ mutagene signature identify -i ../../samples/PD3851a.vcf -g ../../samples/hg38.2bit -f VCF -s10 -b``

---------------
3.4.2.Identify Output
---------------

======  =========  =======  =====  ===========  ===========  ==========  ============
sample  signature  exp       mut    exp_CI_low  exp_CI_high  mut_CI_low  mut_CI_high
======  =========  =======  =====  ===========  ===========  ==========  ============
VCF     2          0.0948    81     0.0921      0.0975       79          83
VCF     3          0.0383    33     0.0340      0.0427       29          36
VCF     4          0.0109    9      0.0081      0.0136       7           12
VCF     5          0.0746    64     0.0722      0.0769       62          66
VCF     6          0.1468    125    0.1417      0.1518       121         129
VCF     7          0.0182    16     0.0154      0.0210       13          18
VCF     8          0.0539    46     0.0505      0.0572       43          49
VCF     9          0.0291    25     0.0264      0.0318       23          27
VCF     10         0.0138    12     0.0112      0.0164       10          14
======  =========  =======  =====  ===========  ===========  ==========  ============

*3.5. Calculate the mutational profile for PD3851a.vcf using -g hg38.2bit and MutaGene-10 signature set. Use the bootstrap to calculate confidence intervals and do not account for unexplained variance (non-context dependent mutational processes and unknown signatures):*
-------
3.5.1.Command
-------

``$ mutagene signature identify -i ../../samples/PD3851a.vcf -g ../../samples/hg38.2bit -f VCF -s10  -U -b``

---------------
3.5.2.Identify Output
---------------

====== ========= =======    =====  ===========  ===========  ========== ============
sample signature exp        mut     exp_CI_low  exp_CI_high  mut_CI_low  mut_CI_high
====== ========= =======    =====  ===========  ===========  ========== ============
VCF     2        0.1233     105     0.1203      0.1263       103         108
VCF     3        0.1987     170     0.1944      0.2030       166         173
VCF     4        0.0697     59      0.0676      0.0717       58          61
VCF     5        0.0878     75      0.0850      0.0906       73          77
VCF     6        0.1820     155     0.1782      0.1858       152         159
VCF     7        0.0980     84      0.0956      0.1005       82          86
VCF     8        0.1047     89      0.1016      0.1077       87          92
VCF     9        0.0633     54      0.0608      0.0658       52          56
VCF     10       0.0708     60      0.0682      0.0735       58          63
====== ========= =======    =====  ===========  ===========  ========== ============

*3.6. Calculate the mutational profile for PD3851a.vcf using -g hg38.2bit and MutaGene-10 signature set and frobenius method:*
--------
3.6.1.Command
-------

``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m frobenius``

---------------
3.6.2.Identify Output
---------------

======  ==========  ==========  ===========
sample  signature    exposure    mutations
======  ==========  ==========  ===========
VCF     5            1.0000      853
======  ==========  ==========  ===========

*3.7. Calculate the mutational profile for PD3851a.vcf using -g hg38.2bit and MutaGene-10 signature set and frobeniuszero method:*
-------
3.7.1.Command
-------

``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m frobeniuszero``

---------------
3.7.2.Identify Output
---------------

======  ==========  ==========  =========
sample  signature   exposure    mutations
======  ==========  ==========  =========
VCF     5           0.7376      629
======  ==========  ==========  =========

3.8. Calculate the mutational profile for PD3851a.vcf using -g hg38.2bit and MutaGene-10 signature set and either the js or divergencejs method:*
-------
3.8.1.Command
-------

``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m js``


``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m divergencejs``

----------------
3.8.2.Identify Output
----------------

Both methods generate the same output for this input

======  =========  ==========  ==========
sample  signature   exposure    mutations
======  =========  ==========  ==========
VCF     2           0.0795     68
VCF     3           0.1634     139
VCF     4           0.0244     21
VCF     5           0.0756     64
VCF     6           0.2012     172
VCF     7           0.0791     67
VCF     8           0.0756     64
VCF     9           0.1020     87
VCF     10          0.1186     101
======  =========  ==========  ==========

3.9. Calculate the mutational profile for PD3851a.vcf using -g hg38.2bit and MutaGene-10 signature set and either the compat, aicc, bic or aiccz method:*
-------
3.9.1.Command
-------

``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m compat``


``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m aic``


``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m bic``


``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m aiccz``

---------------
3.9.2.Identify Output
---------------
All 4 methods generate the same output for this input

======  =========  ==========  =========
sample  signature   exposure   mutations
======  =========  ==========  =========
VCF     2           0.0973     83
VCF     3           0.0536     46
VCF     5           0.0825     70
VCF     6           0.1687     144
VCF     7           0.0220     19
VCF     8           0.0296     25
VCF     9           0.0213     18
VCF     10          0.0034     3
======  =========  ==========  =========
