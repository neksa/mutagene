"""Verify numerical diagnostics and the optional objective scaling."""

import numpy as np

from mutagene.signatures.identify import (
    NegLogLik,
    NegLogLikScaled,
    decompose_mutational_profile_counts,
)


def test_scaled_objective_preserves_the_optimum_definition():
    matrix = np.array([[0.8, 0.2], [0.2, 0.8]])
    counts = np.array([80, 20])
    weights = np.array([0.9, 0.1])
    assert np.isclose(
        NegLogLikScaled(weights, matrix, counts), NegLogLik(weights, matrix, counts) / 100
    )


def test_scaled_method_reports_optimizer_diagnostics():
    matrix = np.full((96, 1), 1 / 96)
    counts = np.full(96, 500)
    diagnostics = {}
    weights, _, _ = decompose_mutational_profile_counts(
        counts,
        (matrix, ["uniform"]),
        func="MLE-scaled",
        enable_dummy=False,
        diagnostics=diagnostics,
    )
    assert diagnostics["optimizer_success"]
    assert diagnostics["fallback"] is None
    assert np.isclose(weights[0], 1)
