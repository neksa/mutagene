# Select sparse signature mixtures

Experimental subset selection penalizes unnecessary known signatures while refitting their likelihoods.
The default `MLE` method remains unchanged.

## Run a sparse fit

AICc requires enough mutations for the model dimensions and its small-sample correction.
BIC also requires adequate mutation support for useful interpretation.
Neither criterion supplies a posterior probability that a signature is present.

```sh
mutagene signature -f MATRIX -i matrix.tsv --method sparse-bic \
  -o sparse.tsv --diagnostics-out sparse.json
```

Use `--method sparse-aicc` to compare AICc.
Use `--no-unexplained-variance` to remove the context-independent background components.
The command retains the existing option name for compatibility.

## Understand the search

The fitter minimizes multinomial negative log likelihood under nonnegative weights that sum to one.
It scales the objective by mutation burden and supplies its analytic gradient to SLSQP.
It refits each candidate subset before comparing information criteria.
Backward steps remove known signatures; forward steps reconsider excluded signatures.
The search stops when no single removal or addition improves the criterion.
This procedure finds a local subset optimum, not a guaranteed global optimum.

The information criteria use unscaled negative log likelihood.
For `k` known signatures and six available background components, the simplex has `k + 5` nominal dimensions.
Without background components, it has `k - 1` nominal dimensions.
Mixture boundaries and correlated signatures limit the usual information-criterion interpretation.
Validation must measure attribution accuracy independently.

The fitter retains all six background components as nuisance parameters during known-signature subset selection.
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

Diagnostics include selected known-signature indices, subset-refit count, criterion value, and the search trace.
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
