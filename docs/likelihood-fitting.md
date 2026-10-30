# Fit and validate signature likelihoods

Python, CLI, and web interfaces retain full-model MLE as the default.
Use `func="MLE"` or `--method MLE` for full-model likelihood fitting.
Signature catalog choices remain separate from fitting-method choices.

## Understand the probability model

The fitting boundary validates finite nonnegative counts and signature entries.
It normalizes each positive signature column without modifying the caller's matrix.
It rejects zero columns, empty profiles, and observed channels outside the reference's support.
Background components can supply support when the known signatures cannot explain an observed channel.

MLE minimizes multinomial negative log likelihood per mutation.
Dividing by mutation burden preserves its optimum.
Weights remain nonnegative and sum to one.
The solver uses analytic derivatives and an explicit equality constraint.
The objective does not conditionally normalize weights during derivative calculations.

`MLE-scaled`, `mlez`, and `Compat` use the same repaired likelihood solver.
`Compat` no longer applies its historical discontinuous objective.
Existing `AICc`, `BIC`, and other non-MLE objectives retain their separate implementations.
Their diagnostic records identify any normalized-NNLS fallback.

## Inspect convergence

Use `--diagnostics-out diagnostics.json` with the CLI.

The solver validates a feasible candidate before accepting it.
`optimality_gap` bounds the remaining negative log likelihood improvement per mutation.
The acceptance tolerance is `2e-5` nats per mutation.
This bound does not quantify uncertainty in signature weights or biological attribution.

`solver_success` records the optimizer's own status.
`optimizer_success` records successful feasibility and optimality validation.
A candidate can pass validation despite an unsuccessful solver status.
Conversely, a successful solver status does not bypass validation.
`candidate_renormalized` records any final feasibility correction.

An uncertified candidate triggers a likelihood retry from a different initialization.
After three unsuccessful attempts, fitting raises `RuntimeError`.
MLE and sparse refits never replace unresolved likelihood fits with an NNLS estimate.
NNLS can supply an initialization; it does not determine the returned MLE estimate.

## Compare likelihoods correctly

Every fitting method reports `LogLik` using the original counts and normalized signature columns.
The value omits the multinomial coefficient, which remains constant for a given count vector.
Compare models for the same count vector using this value.
Do not compare raw `LogLik` values across different mutation burdens as a fit-quality ranking.

## Separate convergence from attribution

Accurate optimization does not guarantee accurate signature attribution.
Correlated signatures and sampling noise can support multiple plausible exposure patterns.
Sparse methods add explicit subset selection to the shared likelihood solver.
Each sparse refit must pass the same feasibility and optimality checks.
The [sparse-fitting guide](sparse-signatures.md) explains the default selection policy.
