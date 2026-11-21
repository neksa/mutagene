Estimate signature exposures
============================

Use ``mutagene signature`` to fit a sample against known signatures.
The command reports an exposure and estimated mutation count for each reported signature.
The original MutaGene study uses nonsmooth nonnegative matrix factorization to derive MutaGene-5 and MutaGene-10.
This command holds the selected signatures fixed and estimates their exposures in new samples.
See the `original study <https://doi.org/10.1093/nar/gkx367>`_.

Fit mutations or sample counts
------------------------------

For mutation calls, choose the reference assembly that matches their coordinates.
For VCF input, select a single tumor sample before analysis.

.. code-block:: bash

   mutagene signature -i sample1.maf -g hg19 -s COSMICv2 -o exposures.tsv

For a labelled SBS96 count matrix, use ``-f MATRIX``.
The bundled demo supplies a small matrix for an offline run.

.. code-block:: bash

   mutagene demo --output-dir demo
   mutagene signature -f MATRIX -i demo/matrix.tsv -o exposures.tsv --diagnostics-out diagnostics.json

Fit pentanucleotide and indel signatures
----------------------------------------

Select ``--context SBS1536`` for pentanucleotide signatures or ``--context ID83`` for small insertions and deletions.
Supply a reference TSV with one signature per column.
The first column must contain channel labels under ``Type`` or ``MutationType``.
The count matrix uses the same labels and one sample per column.
MutaGene aligns labels before fitting. Channel order can differ between files.
Each reference column requires finite nonnegative values with positive total mass.
The fitter normalizes the column to unit mass.

For SBS1536, each channel records two flanking bases on each side, such as ``AA[C>T]GT``.
The CLI accepts variant calls with a matching assembly, or an existing SBS1536 matrix.
To export the per-sample matrix, use the profile command.

.. code-block:: bash

   mutagene profile --per-sample --context SBS1536 -g hg19 -i tumors.maf -o penta-counts.tsv
   mutagene signature --context SBS1536 -f MATRIX -i penta-counts.tsv --reference-matrix penta-signatures.tsv -o exposures.tsv

For ID83, supply a classified count matrix from SigProfilerMatrixGenerator or another compatible source.
The ID83 schema distinguishes indel size, repeat context, and deletion microhomology.
For example, ``2:Del:M:1`` describes a two-base deletion with one base of microhomology.

.. code-block:: bash

   mutagene signature --context ID83 -f MATRIX -i indel-counts.tsv --reference-matrix indel-signatures.tsv -o exposures.tsv

Both extended fits use the supplied catalog without the six SBS96 background components.
The default sparse BIC search, diagnostics, and seeded bootstrap apply to both formats.
The small-indel classifier runs upstream of this command.
The bundled reference collections use SBS96.

Choose a signature collection
-----------------------------

The CLI defaults to ``COSMICv2``.
Use ``--signatures`` or ``-s`` to select another collection.

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

Select a fitting method
-----------------------

The default ``sparse-bic`` selects a subset of reference signatures and background components.
For full-model likelihood fitting, use ``--method MLE``.
The repaired likelihood solver can change estimates relative to older releases.

Use ``--keep-only`` to restrict the reference to named signatures.
Use ``--bootstrap`` to estimate confidence intervals by resampling mutations within each sample.
The default uses 100 replicates and a 95% percentile interval.

.. code-block:: bash

   mutagene signature -f MATRIX -i demo/matrix.tsv -s MGA --bootstrap --seed 42 -o bootstrap.tsv

Common arguments
----------------

``--keep-filtered`` and ``--filter-column`` control which variants are read, and
``--params-out`` and ``--params-in`` record and replay a run. See
:doc:`common_options`.

Labelled count matrices and diagnostics
---------------------------------------

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
Each accepted likelihood fit meets an optimality tolerance of ``2e-5`` nats per mutation.
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
