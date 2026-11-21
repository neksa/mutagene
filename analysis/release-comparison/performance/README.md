# Fixed-catalog attribution accuracy and runtime

The [publication figure](performance.pdf) compares five attribution tools with a nonnegative least squares (NNLS) baseline.
Each point shows one profile. Each asterisk marks the lowest observed mean in its column.
MutaGene has the lowest mean error on the 96-profile panel.
The NNLS baseline has the shortest median batch time on the 24-profile timing subset.
The original timing pilot covered only MutaGene and SigProfilerAssignment on clean two-signature mixtures.

| Method | Mean attribution error, 96 profiles | Median seconds, 24 profiles |
|---|---:|---:|
| MutaGene | 0.1935 | 18.75 |
| SigProfilerAssignment | 0.2840 | 4.30 |
| deconstructSigs | 0.3126 | 61.27 |
| signeR with background | 0.3333 | 871.30 |
| MuSiCal | 0.3736 | 4.05 |
| signeR | 0.4070 | 701.73 |
| NNLS baseline | 0.4529 | 0.70 |

MutaGene has the lowest mean error in nine of twelve burden and background combinations.
deconstructSigs leads at 2,000 and 50,000 mutations with missing reference signatures.
At 50,000 mutations, its mean error is 0.1413 and MutaGene's is 0.2501.
On clean profiles with 2,000 mutations, signeR with background has the lowest mean error (0.0872).
Each combination includes eight profiles.

## What each method fits

Each method receives the same 67 COSMIC v3 signatures and the same 96 labelled SBS channels.

| Method | Fit and background treatment |
|---|---|
| MutaGene | Joint sparse BIC selects known signatures and six substitution-specific background components. |
| SigProfilerAssignment 1.1.5 | The custom-catalog workflow uses default assignment settings and one worker. [Method](https://pmc.ncbi.nlm.nih.gov/articles/PMC10746860/) |
| MuSiCal 1.0.0 | The upstream `refit` function uses its default bidirectional likelihood selection. [Method](https://doi.org/10.1038/s41588-024-01659-0) |
| deconstructSigs 1.9.0 | The upstream R fitter uses its default 6% cutoff and reports the remaining mass as unknown. [Method](https://doi.org/10.1186/s13059-016-0893-4) |
| NNLS baseline | SciPy minimizes squared reconstruction error across all catalog columns with nonnegative weights. [Solver](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html) |
| signeR 2.14.0 | The archived Bayesian fits use a fixed catalog, with and without six background components. [Method](https://doi.org/10.1093/bioinformatics/btw572) |

The benchmark scales NNLS, MuSiCal, and SigProfilerAssignment activities to sum to one.
These methods receive zero aggregate background mass in scoring.
For deconstructSigs, the score counts the unknown fraction as aggregate background mass.
This fraction measures discarded weight and has no specific biological interpretation.
MutaGene and the extended signeR configuration estimate background weights directly.

## Profiles and error

The panel includes eight empirical cancer exposure distributions, three mutation burdens, and four background conditions.
Each combination contributes one synthetic profile from the [matched signeR panel](../signer-panel-errors.tsv).
The source generator uses seed 202610273 and COSMIC v3.3.1.
The fitting catalog uses COSMIC v3.
The conditions use clean profiles, missing reference signatures, uniform background, and substitution-specific background.
Each contaminated profile contains 20% additional generating mass.
Two conditions match MutaGene's background family.
The same profiles appear in earlier benchmark work.
Independent tumor cohorts require separate validation.

[Counts](counts.tsv), [reference probabilities](reference.tsv), [generating weights](truth.tsv), and [sample metadata](samples.tsv) specify the complete experiment.
The [input checksums](../../local-artifacts.json) identify the source files under `signature-benchmark-2026/background-results/fresh`.
We select `sample_index == 0` in each cancer, burden, and condition combination.
Known weights follow fitting-catalog order.
The final truth column contains the total generating mass outside the known signatures.

[Per-profile errors](errors.tsv) give half the summed absolute deviations across known-signature weights and aggregate background mass.
Error ranges from zero for exact recovery to one for disjoint mixtures.
The table below the box plots reports arithmetic means.
Asterisks denote descriptive minima using full precision.
Boxes cover the interquartile range. Their central lines mark medians.
Each burden contributes 32 profiles.
The [archived signeR results](../signer-panel-errors.tsv) supply both signeR rows.
Nineteen of twenty evaluable independent chain pairs exceed 1.01 for maximum active-component R-hat.
The [signeR assessment](../../signature-benchmark-2026/signer-report.qmd) documents convergence and sampler settings.

## Accuracy by background condition

Each cell gives mean attribution error across 24 synthetic profiles and three burdens.

| Method | Clean | Missing signatures | Uniform noise | Substitution-specific noise |
|---|---:|---:|---:|---:|
| MutaGene | 0.091 | 0.318 | 0.185 | 0.181 |
| SigProfilerAssignment | 0.124 | 0.359 | 0.294 | 0.360 |
| deconstructSigs | 0.318 | 0.285 | 0.276 | 0.371 |
| MuSiCal | 0.234 | 0.385 | 0.403 | 0.472 |
| signeR | 0.287 | 0.445 | 0.399 | 0.497 |
| signeR with background | 0.295 | 0.411 | 0.329 | 0.298 |
| NNLS baseline | 0.316 | 0.474 | 0.468 | 0.554 |

The two noise controls add mass within MutaGene's background family.
The missing-signature control allocates 10% to each of two signatures outside the fitting catalog.
Medo et al. also test individual signatures, between-group activity differences, reference pruning, and graded missing-signature fractions.
They assess real-tumor downsampling and agreement between fitters.
Our [PCAWG example](../cohort-examples/README.md) fits the original high-burden catalogs.
A matched downsampling comparison and graded missing-signature sweep remain open. [Published benchmark](https://doi.org/10.1038/s41467-024-53711-6)

## Runtime

The [timing records](times.tsv) include one warm-up and three measured rounds.
Each round uses the same 24 profiles on one Apple M1 Max running macOS 27.2.
Each burden and condition combination contributes two profiles.
Each cancer exposure distribution contributes three profiles.
The timing column in the sample metadata identifies the profiles before fitting.

Each method runs sequentially in a fresh process with one worker and one numerical library thread.
The runner rotates method order between rounds.
Elapsed time includes imports, input parsing, fitting, and native and common output writes.
SigProfilerAssignment disables plots and writes its native assignment files.
The deconstructSigs adapter sources the upstream R fitting functions directly.
MutaGene uses the Python fitting interface for the supplied reference.
These measurements describe batch workflows on this machine.
The earlier 32-profile pilot times the MutaGene CLI and uses different inputs.
The archived signeR jobs use concurrent workers.
Across three measured batches, standard signeR takes 683.61 to 857.45 seconds.
With background, it takes 862.69 to 876.26 seconds.
The separate timing runner uses the same 24-profile matrix, hardware, thread limits, and repetition count as the other methods.
It starts one fresh R process per batch and alternates the two configurations.
Each sample uses its archived primary-chain seed, 10,000 burn-in iterations, and 2,000 retained draws.
Timing includes package startup, fitting, median exposure extraction, and writing the weight matrix.
Each timed round repeats the same seeds. The repeats measure timing variability.
The archived independent chains provide the convergence diagnostics.
The [repeat-fit errors](signer-timing-errors.tsv) record the fresh estimates on the 24 timing profiles.
Fresh signeR scores on the timing subset differ from the archived estimates.
Mean error changes from 0.4042 to 0.4181 for the standard fit.
With background, it changes from 0.3372 to 0.3288.
Individual absolute error differences reach 0.328 and 0.323, respectively.
The figure uses the complete 96-profile accuracy panel for every method.
Median batch time falls from 19.48 seconds to 18.75 seconds, a 3.7% reduction.
NNLS takes 0.70 seconds, MuSiCal takes 4.05 seconds, and SigProfilerAssignment takes 4.30 seconds.
deconstructSigs takes 61.27 seconds.

## Why MutaGene takes time

Joint sparse BIC compares many neighboring subsets.
Each feasible subset can require a constrained likelihood optimization with SLSQP.
On the 32-profile pilot, the previous fitter calls the likelihood fitter 8,029 times and SLSQP 6,977 times.
In the instrumented run, likelihood fitting accounts for 6.89 of 8.89 seconds.
Profiling adds overhead.

For each proposed addition, the optimized fitter bounds the best possible improvement in log likelihood.
The tangent at the current mixture bounds the concave log likelihood across the enlarged simplex.
A linear function reaches its maximum at a simplex vertex.
If even that bound cannot offset the BIC penalty, the fitter skips the refit.
The search retains borderline candidates and checks the remaining fits under the existing optimality tolerance.
The [pilot diagnostics](pilot-refits.tsv) give 5,550 refits after this change.
Both versions select identical subsets and produce identical exposure tables on these 32 profiles.
They also produce byte-identical known-signature weights and aggregate background estimates on the 96-profile panel.
After optimization, the 16 profiles with 100 mutations require 13,200 solver iterations.
The 16 profiles with 2,000 mutations require 32,390 iterations.
The larger constrained fits provide the next profiling target.

The [repeated CLI pilot](pilot-times.tsv) uses the original 32 clean profiles.
Median MutaGene time falls from 5.83 to 5.32 seconds, a 8.8% reduction.
SigProfilerAssignment takes 3.53 seconds on the same inputs.
Each workflow uses three measured fresh processes after one warm-up.

## Reproduce the comparison

Use Python 3.12.12, R 4.6.2, Git, and uv.
Run these commands from the repository root.
The environment installs the previous MutaGene implementation.
The runner imports the optimized implementation from this checkout.

```sh
git worktree add --detach /tmp/mutagene-before d9eac3fe4c28c1c36cfcaecd5fba47ffcf252303
git clone https://github.com/raerose01/deconstructSigs /tmp/deconstructSigs
git -C /tmp/deconstructSigs checkout 41a705c5d80848121347d448cf9e2c58ff6b81ac
uv venv --python 3.12.12 /tmp/mutagene-performance
uv pip install --python /tmp/mutagene-performance/bin/python /tmp/mutagene-before \
  -r analysis/signature-benchmark-2026/requirements-benchmark.txt \
  matplotlib==3.11.2 seaborn==0.13.2 \
  'musical @ git+https://github.com/parklab/MuSiCal@2b0c246972376911893b65868c1b756671ef3436'
cd analysis/release-comparison/performance
/tmp/mutagene-performance/bin/python benchmark.py timing /tmp/deconstructSigs /tmp/attribution-timing
/tmp/mutagene-performance/bin/python benchmark.py all /tmp/deconstructSigs /tmp/attribution-errors
/tmp/mutagene-performance/bin/python figure.py
```

The last command renders the figure from the published tables in this directory.
To render a new run, pass its error table with `--errors` and its output filename with `--output`.
The runner checks sample identity, signature identity, unit mass, repeated weights, and agreement between MutaGene versions.
For the CLI pilot, generate the inputs with [compare.py](../compare.py).
Pass its output directory and a new results directory to `time_pilot.py`.

The [preprint source](../../../manuscript/preprint.tex) includes the vector figure.
From the manuscript directory, run `pdflatex -halt-on-error preprint.tex` twice to build the draft.

### Reproduce the signeR timing

Install signeR 2.14.0 from Bioconductor 3.23 in an isolated R library.
The timing run uses the macOS arm64 R 4.6.2 framework with its reference BLAS.
The runner checks the installed signeR version.

From this directory, run the following command with the isolated library on the R library path.

```sh
R_LIBS=/path/to/library /tmp/mutagene-performance/bin/python time_signer.py \
  /tmp/attribution-signer --rscript /usr/local/bin/Rscript
```

The [seeds](signer-seeds.tsv) select the archived primary fit for each configuration.
The [background catalog](signer-background.tsv) contains the same six substitution-specific components as the accuracy study.
