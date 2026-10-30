"""Fit and certify multinomial likelihoods on the probability simplex."""

import numpy as np
from scipy.optimize import minimize, nnls

DEFAULT_METHOD = "sparse-bic"
DEFAULT_BACKGROUND_SELECTION = "joint"
OPTIMALITY_TOLERANCE = 2e-5


def validate_inputs(matrix, counts):
    """Validate counts and convert each signature column to a probability distribution."""
    matrix = np.asarray(matrix, dtype=float)
    counts = np.asarray(counts, dtype=float)
    if matrix.ndim != 2 or counts.ndim != 1 or matrix.shape[0] != len(counts):
        raise ValueError("Signature matrix and count dimensions differ")
    if matrix.shape[1] == 0:
        raise ValueError("Provide at least one signature")
    if (
        not np.isfinite(matrix).all()
        or not np.isfinite(counts).all()
        or (matrix < 0).any()
        or (counts < 0).any()
    ):
        raise ValueError("Provide finite nonnegative signatures and counts")
    totals = matrix.sum(axis=0)
    burden = counts.sum()
    if not np.isfinite(totals).all() or (totals <= 0).any():
        raise ValueError("Each signature must have positive finite total probability")
    if not np.isfinite(burden) or burden <= 0:
        raise ValueError("Profile is all zeros or has an invalid total count")
    return matrix / totals, counts


def negative_log_likelihood(weights, matrix, counts):
    """Return count NLL without the constant multinomial coefficient or implicit weight normalization."""
    counts = np.asarray(counts, dtype=float)
    observed = counts > 0
    predicted = np.asarray(matrix, dtype=float)[observed] @ np.asarray(weights, dtype=float)
    if not np.isfinite(predicted).all() or (predicted <= 0).any():
        return np.inf
    return float(-np.dot(counts[observed], np.log(predicted)))


def likelihood_objective(weights, design, target):
    """Return per-mutation NLL and its analytic gradient on observed channels."""
    predicted = design @ weights
    # Invalid trial points have infinite NLL; only the derivative uses a finite denominator.
    denominator = np.where(predicted > 0, predicted, 1e-300)
    with np.errstate(over="ignore", invalid="ignore"):
        gradient = -((design / denominator[:, None]).T @ target)
    if (predicted <= 0).any() or not np.isfinite(predicted).all():
        return np.inf, gradient
    return float(-np.dot(target, np.log(predicted))), gradient


def fit_likelihood(matrix, counts, diagnostics=None):
    """Return a feasible likelihood fit with a certified normalized objective gap.

    The gap bounds remaining NLL improvement per mutation, not uncertainty in signature weights.
    A solver status alone does not determine acceptance.
    """
    matrix, counts = validate_inputs(matrix, counts)
    observed = counts > 0
    design = matrix[observed]
    if (design.sum(axis=1) == 0).any():
        raise ValueError("The signature reference cannot explain an observed channel")
    target = counts[observed] / counts.sum()
    size = matrix.shape[1]
    uniform = np.full(size, 1 / size)
    try:
        start, _ = nnls(matrix, counts / counts.sum())
        start = start / start.sum() if start.sum() else uniform.copy()
    except RuntimeError:
        start = uniform.copy()
    start = start * (1 - 1e-6) + 1e-6 * uniform
    attempts = []
    best = uniform.copy()
    best_nll = np.inf
    for attempt in range(3):
        if attempt == 1:
            start = uniform.copy()
        elif attempt == 2:
            start = best * (1 - 1e-8) + 1e-8 * uniform
        if size == 1:
            raw, solver_success, message, iterations = uniform, True, "Single signature", 0
        else:
            result = minimize(
                likelihood_objective,
                start,
                args=(design, target),
                jac=True,
                method="SLSQP",
                bounds=[(0, 1)] * size,
                constraints={
                    "type": "eq",
                    "fun": lambda w: w.sum() - 1,
                    "jac": lambda w: np.ones(len(w)),
                },
                options={"ftol": 1e-12 if attempt == 0 else 1e-14, "maxiter": 2000},
            )
            raw = np.asarray(result.x, dtype=float)
            solver_success, message, iterations = (
                bool(result.success),
                str(result.message),
                int(result.nit),
            )
        gap, nll, weights = np.inf, np.inf, None
        if raw.shape == (size,) and np.isfinite(raw).all() and raw.min() >= -1e-8 and raw.sum() > 0:
            weights = np.maximum(raw, 0)
            weights /= weights.sum()
            nll, gradient = likelihood_objective(weights, design, target)
            if np.isfinite(nll) and np.isfinite(gradient).all():
                gap = max(0.0, float(gradient @ weights - gradient.min()))
                if nll < best_nll:
                    best, best_nll = weights, nll
        accepted = bool(np.isfinite(nll) and gap <= OPTIMALITY_TOLERANCE)
        attempts.append(
            {
                "solver_success": solver_success,
                "message": message,
                "iterations": iterations,
                "accepted": accepted,
            }
        )
        if accepted:
            if diagnostics is not None:
                diagnostics.update(
                    optimizer_success=True,
                    optimizer_message=message,
                    solver_success=solver_success,
                    optimizer_iterations=sum(item["iterations"] for item in attempts),
                    optimizer_retries=attempt,
                    optimizer_attempts=attempts,
                    optimality_gap=gap,
                    optimality_tolerance=OPTIMALITY_TOLERANCE,
                    nll_per_mutation=nll,
                    weight_sum=float(weights.sum()),
                    estimator="multinomial_mle",
                    fallback=None,
                    candidate_renormalized=bool(np.any(raw < 0) or abs(raw.sum() - 1) > 1e-12),
                )
            return weights
    if diagnostics is not None:
        diagnostics.update(
            optimizer_success=False,
            optimizer_attempts=attempts,
            optimizer_retries=2,
            estimator="multinomial_mle",
            fallback=None,
        )
    raise RuntimeError(
        f"Likelihood optimization fails certification after three attempts: {message}"
    )
