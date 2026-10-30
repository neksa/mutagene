Estimate signature exposures
============================

``mutagene signature`` estimates mutational signature exposures from mutation data.
If you install MutaGene in a virtual environment, activate that environment before running commands.

Inspect command arguments:

.. code-block:: bash

   mutagene signature identify -h

2. Arguments
------------

**2.1.Command:** mutagene signature identify [arguments]

**2.2.Required Arguments (must be specified):**


Choose one of five signature collections with ``--signatures`` (``-s``):

==================  ==========  ==============================================
Value               Alias       Collection
==================  ==========  ==============================================
MGA                 5           MutaGene-5, 5 signatures derived by MutaGene
MGB                 10          MutaGene-10, 10 signatures derived by MutaGene
COSMICv2            30          COSMIC v2, 30 signatures
COSMICv3            49          COSMIC v3, 67 SBS signatures (SBS1-SBS85)
KUCAB               53          54 profiles of environmental mutagen exposure,
                                including a Control profile
==================  ==========  ==============================================

Numeric aliases are historical names, not signature counts.
``-s 49`` loads 67 COSMIC v3 signatures.
``-s 53`` loads 54 KUCAB profiles.

Read more about the MutaGene signature package at https://www.ncbi.nlm.nih.gov/research/mutagene/signatures#mutational_signatures.
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
                            (default method is sparse-bic)\ :sup:`1`
-m [METHOD]                 Short form of --method [METHOD] argument                             -m frobenius
--no-unexplained-variance   not to account for unexplained variance from context-independent     --no-unexplained-variance
                            mutational processes and unknown signatures
-U                          Short form of --no-unexplained-variance argument                     -U
--bootstrap                 Use the bootstrap to calculate 95% confidence intervals              --bootstrap
                            for VCF files
-b                          Short form of --bootstrap argument                                   -b
==========================  ===================================================================  =========================

Bootstrap sampling draws mutations with replacement from each sample 100 times.
The output reports mean estimates and 95% confidence intervals.

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
compat             NegLogLik (compatibility alias for repaired MLE)
aicc               AICc (AIC corrected for small samples)
bic                BIC
aiccz              BIC (BIC with added context-independent signatures)
================= =====================================================================================


3. Examples
-----------
Fit MutaGene-10 signatures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fit MutaGene-10 signatures in ``PD3851a.vcf`` with hg38 and ``mlez``.
The ``mle`` method gives the same output for this input.

3.1.1.Command
-------------

``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 --method mlez``

3.1.2.Identify Output
---------------------

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


Fit MutaGene-5 signatures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fit MutaGene-5 signatures in ``sample1.maf`` with ``hg38.2bit``.

3.2.1.Command
-------------

``$ mutagene signature identify -i sample1.maf -g hg38.2bit -s5``

3.2.2.Identify Output
---------------------

=============================  ===========  ==========  ===========
sample                          signature    exposure    mutations
=============================  ===========  ==========  ===========
TCGA-50-6593-01A-11D-1753-08    2            0.0348      5
TCGA-50-6593-01A-11D-1753-08    3            0.0691      11
=============================  ===========  ==========  ===========

Save exposures to a file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fit MutaGene-10 signatures and save exposures to ``out.tsv``.
This example reads input and genome files from the samples directory.

3.3.1.Command
-------------

``$ mutagene signature identify -i ../../samples/PD3851a.vcf -g ../../samples/hg38.2bit -f VCF -s10 -o ../../out/out.tsv``

3.3.2.Identify Output
---------------------
As for example (3.1) except sent to file out.tsv instead of screen.

Estimate confidence intervals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use bootstrap sampling to estimate confidence intervals for MutaGene-10 exposures.

3.4.1.Command
-------------

``$ mutagene signature identify -i ../../samples/PD3851a.vcf -g ../../samples/hg38.2bit -f VCF -s10 -b``

3.4.2.Identify Output
---------------------

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

Bootstrap without unexplained variance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Estimate bootstrap intervals without the unexplained-variance components.

3.5.1.Command
-------------

``$ mutagene signature identify -i ../../samples/PD3851a.vcf -g ../../samples/hg38.2bit -f VCF -s10  -U -b``

3.5.2.Identify Output
---------------------

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

Use Frobenius fitting
~~~~~~~~~~~~~~~~~~~~~~~~~~

Fit MutaGene-10 signatures with the Frobenius method.

3.6.1.Command
-------------

``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m frobenius``

3.6.2.Identify Output
---------------------

======  ==========  ==========  ===========
sample  signature    exposure    mutations
======  ==========  ==========  ===========
VCF     5            1.0000      853
======  ==========  ==========  ===========

Use FrobeniusZero fitting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fit MutaGene-10 signatures with the FrobeniusZero method.

3.7.1.Command
-------------

``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m frobeniuszero``

3.7.2.Identify Output
---------------------

======  ==========  ==========  =========
sample  signature   exposure    mutations
======  ==========  ==========  =========
VCF     5           0.7376      629
======  ==========  ==========  =========

Use Jensen–Shannon fitting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fit MutaGene-10 signatures with ``js`` or ``divergencejs``.

3.8.1.Command
-------------

``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m js``


``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m divergencejs``

3.8.2.Identify Output
---------------------

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

Compare legacy methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compare the legacy fitting methods below.

3.9.1.Command
-------------

``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m compat``


``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m aic``


``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m bic``


``$ mutagene signature identify -i PD3851a.vcf -g hg38.2bit -f VCF -s10 -m aiccz``

3.9.2.Identify Output
---------------------
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

Common arguments
----------------

``--keep-filtered`` and ``--filter-column`` control which variants are read, and
``--params-out`` and ``--params-in`` record and replay a run. See
:doc:`common_options`.

Labelled count matrices and diagnostics
---------------------------------------

Use ``-f MATRIX`` for an SBS96 count matrix with ``Type`` and sample columns.
The reader validates channel labels and preserves sample identity.
Matrix input does not require a genome file.

.. code-block:: bash

   mutagene signature -f MATRIX -i matrix.tsv -o signatures.tsv --diagnostics-out diagnostics.json

Diagnostics include optimizer status, likelihood optimality bounds, and six null-component weights.
MLE never replaces a failed likelihood fit with normalized NNLS.
Null components represent context-independent substitution patterns, not arbitrary unknown processes.

Use ``--seed`` with bootstrap sampling for reproducible draws.
The percentile method remains the default.
The approximate t-based method uses bootstrap standard deviation rather than replicate SEM.

Python, CLI, and web interfaces default to ``sparse-bic`` with ``joint`` background selection.
Request ``--method MLE`` or ``func="MLE"`` for full-model likelihood fitting.
``MLE-scaled`` and ``Compat`` use the same repaired likelihood solver.
The solver normalizes signature columns and validates finite nonnegative inputs.
It uses per-mutation likelihood, analytic derivatives, and an explicit simplex constraint.
Each accepted fit has an optimality bound below ``2e-5`` nats per mutation.
After three unsuccessful certification attempts, fitting raises ``RuntimeError``.
Every method reports ``LogLik`` from original counts, omitting the constant multinomial coefficient.

Sparse subset selection
-----------------------

``--method sparse-bic`` and ``--method sparse-aicc`` search known-signature subsets.
Each candidate subset receives a constrained multinomial likelihood refit.
Backward removals and forward additions stop at a local information-criterion optimum.
Joint-background BIC starts from the best feasible single component through 1,000 mutations.
It compares all one-component additions and removals at each step.
Above 1,000 mutations, it retains full-model initialization and the preceding search.
Fixed-background BIC, block selection, and AICc retain full-model initialization.
The diagnostics field ``search_start`` records the initialization.
The default method uses BIC to select known and background components.

The six null components form a context-independent background model.
Each component covers one substitution type uniformly across sixteen flanking contexts.
Use ``--background-selection joint`` to select background components individually.
Use ``--background-selection block`` to compare complete models with and without background.
Use ``--background-selection fixed`` to retain all six components as nuisance parameters.
Joint selection counts retained known and background components minus one for the simplex constraint.
Their fitted weights do not constitute a Bayesian prior or prove technical noise.

Use ``--no-unexplained-variance`` to compare a fit without these components.
Inspect ``--diagnostics-out`` for the selected subset, search trace, and optimizer status.
