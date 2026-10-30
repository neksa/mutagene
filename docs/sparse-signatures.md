# Select sparse signature mixtures

Request `--method sparse-bic` to select signature subsets.
The default background selection policy is `joint`.
Each candidate subset receives a constrained multinomial MLE refit.
Use `--method MLE` for full-model fitting without subset selection.

## Run a sparse fit

AICc requires enough mutations for the model dimensions and its small-sample correction.
BIC also requires adequate mutation support for useful interpretation.
Neither criterion supplies a posterior probability that a signature is present.

```sh
mutagene signature -f MATRIX -i matrix.tsv --method sparse-bic \
  -o sparse.tsv --diagnostics-out sparse.json
```

Use `--method sparse-aicc` to compare AICc.
Subset searches cost more computation than a single full-model MLE fit.
Bootstrap sampling repeats subset selection for every replicate.
Use `--no-unexplained-variance` to remove the context-independent background components.
The command retains the existing option name for compatibility.

## Understand the search

The fitter minimizes multinomial negative log likelihood under nonnegative weights that sum to one.
It scales the objective by mutation burden and supplies its analytic gradient to SLSQP.
It refits each candidate subset before comparing information criteria.
Backward steps remove optional components; forward steps reconsider excluded components.
The search stops when no single removal or addition improves the criterion.
This procedure finds a local subset optimum, not a guaranteed global optimum.

The information criteria use unscaled negative log likelihood.
Count scaling preserves full-model MLE weights.
It can change BIC selection because likelihood evidence grows faster than the logarithmic complexity penalty.
For `k` known signatures and six available background components, the simplex has `k + 5` nominal dimensions.
Without background components, it has `k - 1` nominal dimensions.
Mixture boundaries and correlated signatures limit the usual information-criterion interpretation.
Validation must measure attribution accuracy independently.

## Fit profiles with few mutations

Joint-background BIC starts from the best feasible single component through 1,000 mutations.
It compares all one-component additions and removals and accepts the largest BIC improvement.
Above 1,000 mutations, it retains full-model initialization and the preceding search.
Fixed-background BIC, block selection, and AICc retain full-model initialization at every burden.
If no single component covers all observed channels, joint BIC also uses full-model initialization.

This threshold follows synthetic tuning, not a biological boundary.
Inspect `search_start` in the diagnostics to identify the selected initialization.
Direct `fit_sparse` calls accept `search_start="full"` to retain the preceding search at any burden.

## Select background components

Use `--background-selection joint` to select background components individually alongside known signatures.
Use `--background-selection block` to compare complete models with and without all six background components.
Use `--background-selection fixed` to retain all six background components during known-signature selection.
The Python API accepts `background_selection="joint"`, `"block"`, or `"fixed"`.
These policies affect only sparse methods.
`--no-unexplained-variance` removes background components before selection, regardless of policy.

Joint selection counts `k + b - 1` dimensions for `k` retained known signatures and `b` retained background components.
Fixed selection counts all six background components, including components with zero fitted weight.
Block selection compares complete BIC values, including each model's background dimensions.
It excludes the no-background candidate when that reference cannot explain an observed channel.
Numerical failures still stop selection.

Each background component carries a complexity cost only when the candidate model contains that component.
It can select zero known signatures when the background explains the profile.
Such a fit produces no known-signature rows; inspect the background weights in the diagnostics JSON.
It does not apply a Bayesian prior to background mass.

## Interpret background weights

Each background component covers one substitution type uniformly across sixteen flanking contexts.
Their fitted mixture can represent substitution-specific context-independent signal.
Background weight can absorb sampling error, model mismatch, or genuine context-independent mutagenesis.
It does not prove technical noise or quantify arbitrary unknown biology.

The implementation preserves `null_weight`, `null_components`, and existing dummy-signature identifiers for compatibility.
Reports describe these quantities as context-independent background.
A future shrinkage prior would constitute a distinct statistical model.

## Inspect diagnostics

Diagnostics include selected known-signature and background indices, model dimensions, criterion value, and the search trace.
Background indices follow the existing order of the six dummy signatures.
Block diagnostics include both complete-model criterion values when both models have adequate support.
Each accepted step decreases the criterion.
Each refit uses the shared certified likelihood solver.
An uncertified candidate triggers retries from uniform weights and the best feasible candidate.
After three unsuccessful attempts, fitting stops instead of substituting NNLS.
Diagnostics include the maximum optimality gap across all refits.
For sparse fits, retry counts sum across all subset refits.

## Evaluate the result

Compare attribution error, precision, recall, active-signature count, reconstruction, and background recovery together.
A sparse solution alone does not establish correct attribution.
Keep tuning profiles separate from final evaluation profiles.
The benchmark records selection rules and source hashes before evaluation.
