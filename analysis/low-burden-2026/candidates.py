"""Explore low-burden estimators without changing production fitting defaults."""

from math import lgamma

import numpy as np
from scipy.optimize import nnls

from mutagene.signatures import sparse
from mutagene.signatures.likelihood import fit_likelihood, validate_inputs


def fit_ebic(matrix, counts, known_count, gamma):
    """Add a model-space penalty during the existing joint subset search."""
    original = sparse.information_criterion
    size = matrix.shape[1]

    def criterion(nll, dimensions, burden, name):
        selected = dimensions + 1
        log_models = lgamma(size + 1) - lgamma(selected + 1) - lgamma(size - selected + 1)
        return original(nll, dimensions, burden, name) + 2 * gamma * log_models

    # Each benchmark worker runs one candidate at a time in its own process.
    sparse.information_criterion = criterion
    try:
        return sparse.fit_sparse(matrix, counts, known_count, background_selection="joint")
    finally:
        sparse.information_criterion = original


def fit_sparse_nnls(matrix, counts, tolerance, refit=False):
    """Select components using changes in relative Euclidean reconstruction error."""
    matrix, counts = validate_inputs(matrix, counts)
    target = counts / counts.sum()
    norm = np.linalg.norm(target)
    cache = {}

    def fit(selected):
        selected = tuple(sorted(selected))
        if not selected:
            return None
        if selected not in cache:
            weights, residual = nnls(matrix[:, selected], target)
            expanded = np.zeros(matrix.shape[1])
            expanded[list(selected)] = weights
            cache[selected] = (expanded, residual / norm)
        return cache[selected]

    selected = tuple(range(matrix.shape[1]))
    current = fit(selected)
    selected = tuple(np.flatnonzero(current[0] > 1e-7))
    current = fit(selected)
    while True:
        candidates = [
            (
                tuple(i for i in selected if i != remove),
                fit(tuple(i for i in selected if i != remove)),
            )
            for remove in selected
            if len(selected) > 1
        ]
        if candidates:
            subset, candidate = min(candidates, key=lambda item: item[1][1])
            if candidate[1] - current[1] <= tolerance:
                selected, current = subset, candidate
                continue
        additions = [
            (tuple(sorted((*selected, add))), fit((*selected, add)))
            for add in range(matrix.shape[1])
            if add not in selected
        ]
        if not additions:
            break
        subset, candidate = min(additions, key=lambda item: item[1][1])
        if current[1] - candidate[1] <= 2 * tolerance:
            break
        selected, current = subset, candidate
    weights = current[0]
    if refit:
        weights = np.zeros(matrix.shape[1])
        weights[list(selected)] = fit_likelihood(matrix[:, selected], counts)
    else:
        weights = weights / weights.sum()
    return weights


CONFIGS = [
    ("joint-BIC", "bic", 0.0, True, False),
    ("joint-AICc", "aicc", 0.0, True, False),
    *[(f"EBIC-{g}", "ebic", g, True, False) for g in (0.5, 1.0)],
    *[
        (f"NNLS-{t}-{bg}", "nnls", t, bg, False)
        for t in (0.01, 0.02, 0.05, 0.1)
        for bg in (False, True)
    ],
    *[(f"NNLS-MLE-{t}", "nnls", t, True, True) for t in (0.01, 0.02, 0.05)],
]


def fit_candidate(matrix, counts, known_count, config):
    _, method, parameter, background, refit = config
    design = matrix if background else matrix[:, :known_count]
    if method == "forward":
        weights = fit_forward(design, counts, parameter)
    elif method == "nnls":
        weights = fit_sparse_nnls(design, counts, parameter, refit)
    elif method == "ebic":
        weights = fit_ebic(design, counts, known_count, parameter)
    else:
        weights = sparse.fit_sparse(design, counts, known_count, criterion=method)
    return np.pad(weights, (0, matrix.shape[1] - len(weights)))


def fit_forward(matrix, counts, gamma):
    """Start from the best single component and alternate additions with removals."""
    matrix, counts = validate_inputs(matrix, counts)
    size = matrix.shape[1]
    burden = counts.sum()
    cache = {}

    def fit(selected):
        selected = tuple(sorted(selected))
        if not selected or (matrix[counts > 0][:, selected].sum(axis=1) == 0).any():
            return None
        if selected not in cache:
            info = {}
            weights = fit_likelihood(matrix[:, selected], counts, info)
            k = len(selected)
            penalty = (k - 1) * np.log(burden)
            penalty += 2 * gamma * (lgamma(size + 1) - lgamma(k + 1) - lgamma(size - k + 1))
            expanded = np.zeros(size)
            expanded[list(selected)] = weights
            cache[selected] = (expanded, 2 * burden * info["nll_per_mutation"] + penalty)
        return cache[selected]

    starts = [((i,), fit((i,))) for i in range(size)]
    starts = [item for item in starts if item[1] is not None]
    if not starts:
        raise ValueError("No single component explains every observed channel")
    selected, current = min(starts, key=lambda item: item[1][1])
    while True:
        neighbors = [tuple(sorted((*selected, i))) for i in range(size) if i not in selected]
        neighbors += [tuple(i for i in selected if i != remove) for remove in selected]
        candidates = [(s, fit(s)) for s in neighbors]
        candidates = [item for item in candidates if item[1] is not None]
        subset, best = min(candidates, key=lambda item: item[1][1])
        if best[1] >= current[1] - 1e-7:
            return current[0]
        selected, current = subset, best


ROUND2 = [
    ("joint-BIC", "bic", 0.0, True, False),
    *[(f"EBIC-{g}", "ebic", g, True, False) for g in (1.0, 2.0, 4.0)],
    *[(f"forward-EBIC-{g}", "forward", g, True, False) for g in (0.0, 0.5, 1.0, 2.0)],
]
