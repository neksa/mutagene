"""Regressions for the likelihood audit's numerical and probability defects."""

import argparse
import inspect
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from mutagene.cli.signature_menu import SignatureMenu
from mutagene.io.profile import read_signatures
from mutagene.signatures import get_dummy_signatures_lists
from mutagene.signatures.identify import NegLogLik, decompose_mutational_profile_counts
from mutagene.signatures.likelihood import (
    DEFAULT_METHOD,
    OPTIMALITY_TOLERANCE,
    fit_likelihood,
    likelihood_objective,
    validate_inputs,
)

MATRIX = np.array([[0.8, 0.2], [0.2, 0.8]])
COUNTS = np.array([65, 35])


def test_binomial_closed_form_and_analytic_gradient():
    weights = fit_likelihood(MATRIX, COUNTS)
    assert weights == pytest.approx([0.75, 0.25], abs=1e-6)
    probe = np.array([0.6, 0.4])
    _, gradient = likelihood_objective(probe, MATRIX, COUNTS / COUNTS.sum())
    step = 1e-6
    finite = [
        (
            NegLogLik(probe + step * direction, MATRIX, COUNTS)
            - NegLogLik(probe - step * direction, MATRIX, COUNTS)
        )
        / (2 * step * COUNTS.sum())
        for direction in np.eye(2)
    ]
    assert gradient == pytest.approx(finite, abs=1e-8)


def test_objective_has_no_conditional_normalization_kink():
    weights = np.array([0.75, 0.25])
    step = np.array([1e-7, 0])
    center = NegLogLik(weights, MATRIX, COUNTS)
    left = (center - NegLogLik(weights - step, MATRIX, COUNTS)) / step[0]
    right = (NegLogLik(weights + step, MATRIX, COUNTS) - center) / step[0]
    assert left == pytest.approx(-100, abs=1e-4)
    assert right == pytest.approx(left, abs=1e-4)


def test_repeated_counts_preserve_fitted_weights_and_optimality():
    fixture = json.loads((Path(__file__).parent / "data/mle-count-scaling.json").read_text())
    counts = np.array(fixture["counts"])
    matrix, names = read_signatures("COSMICv3")
    fitted = []
    for multiplier in (1, 20, 500):
        info = {}
        weights, _, _ = decompose_mutational_profile_counts(
            counts * multiplier, (matrix, names), diagnostics=info
        )
        assert info["fallback"] is None
        assert info["optimality_gap"] <= OPTIMALITY_TOLERANCE
        assert weights.sum() == pytest.approx(1, abs=1e-12)
        fitted.append(weights)
    assert np.abs(fitted[0] - fitted[2]).sum() < 1e-5
    assert np.abs(fitted[0] - fitted[1]).sum() < 1e-5


@pytest.mark.parametrize("method", ["MLE", "MLE-scaled", "mlez", "Compat"])
def test_likelihood_aliases_use_certified_mle(method):
    info = {}
    weights, _, _ = decompose_mutational_profile_counts(
        COUNTS, (MATRIX, ["A", "B"]), func=method, enable_dummy=False, diagnostics=info
    )
    assert weights == pytest.approx([0.75, 0.25], abs=1e-6)
    assert info["estimator"] == "multinomial_mle"
    assert info["fallback"] is None


@pytest.mark.parametrize(
    "method", ["MLE", "MLE-scaled", "Frobenius", "Cos", "Compat", "sparse-bic"]
)
def test_reported_likelihood_always_uses_original_counts(method):
    weights, _, results = decompose_mutational_profile_counts(
        COUNTS, (MATRIX, ["A", "B"]), func=method, enable_dummy=False
    )
    reported = next(row["score"] for row in results if row["name"] == "LogLik")
    expected = np.dot(COUNTS, np.log(MATRIX @ weights))
    assert reported == pytest.approx(expected, abs=1e-8)


@pytest.mark.parametrize("success", [False, True])
def test_feasible_certified_candidate_overrides_solver_status(monkeypatch, success):
    # The original fallback discards near-optimal candidates with an unsuccessful status.
    monkeypatch.setattr(
        "mutagene.signatures.likelihood.minimize",
        lambda *a, **k: SimpleNamespace(
            x=np.array([0.75, 0.25]) * 1.009, success=success, message="candidate", nit=5
        ),
    )
    info = {}
    weights = fit_likelihood(MATRIX, COUNTS, info)
    assert weights == pytest.approx([0.75, 0.25], abs=1e-12)
    assert info["optimizer_success"]
    assert info["solver_success"] == success
    assert info["candidate_renormalized"]
    assert info["fallback"] is None


def test_successful_but_suboptimal_candidate_requires_retry(monkeypatch):
    candidates = iter([[0.5, 0.5], [0.75, 0.25]])
    monkeypatch.setattr(
        "mutagene.signatures.likelihood.minimize",
        lambda *a, **k: SimpleNamespace(
            x=np.array(next(candidates)), success=True, message="success", nit=1
        ),
    )
    info = {}
    assert fit_likelihood(MATRIX, COUNTS, info) == pytest.approx([0.75, 0.25])
    assert info["optimizer_retries"] == 1
    assert not info["optimizer_attempts"][0]["accepted"]


@pytest.mark.parametrize("candidate", [[0.5, 0.5], [np.nan, 0.5], [-0.1, 1.1], [0, 0]])
def test_unresolved_candidates_raise_without_nnls_fallback(monkeypatch, candidate):
    calls = []

    def broken(*args, **kwargs):
        calls.append(1)
        return SimpleNamespace(x=np.array(candidate), success=True, message="false success", nit=1)

    monkeypatch.setattr("mutagene.signatures.likelihood.minimize", broken)
    with pytest.raises(RuntimeError, match="certification"):
        fit_likelihood(MATRIX, COUNTS)
    assert len(calls) == 3


def test_nnls_initialization_failure_still_solves_likelihood(monkeypatch):
    def broken(*args, **kwargs):
        raise RuntimeError("NNLS initialization fails")

    monkeypatch.setattr("mutagene.signatures.likelihood.nnls", broken)
    assert fit_likelihood(MATRIX, COUNTS) == pytest.approx([0.75, 0.25], abs=1e-6)


def test_column_scaling_does_not_change_signature_probabilities_or_modify_input():
    raw = MATRIX * np.array([0.998577814, 1.000526])
    original = raw.copy()
    info = {}
    weights, _, results = decompose_mutational_profile_counts(
        COUNTS, (raw, ["A", "B"]), enable_dummy=False, diagnostics=info
    )
    assert np.array_equal(raw, original)
    assert weights == pytest.approx([0.75, 0.25], abs=1e-6)
    assert info["reference_column_sum_min"] == pytest.approx(0.998577814)
    ll = next(row["score"] for row in results if row["name"] == "LogLik")
    assert ll == pytest.approx(np.dot(COUNTS, np.log(MATRIX @ weights)))


@pytest.mark.parametrize(
    "matrix,counts",
    [
        (MATRIX, [-1, 101]),
        (MATRIX, [np.nan, 1]),
        (MATRIX, [np.inf, 1]),
        (MATRIX, [0, 0]),
        (MATRIX, [[1, 2]]),
        (MATRIX, [1, 2, 3]),
        ([[1, -1], [0, 2]], COUNTS),
        ([[1, np.nan], [0, 1]], COUNTS),
        ([[1, np.inf], [0, 1]], COUNTS),
        ([[0, 1], [0, 1]], COUNTS),
        (np.empty((2, 0)), COUNTS),
    ],
)
def test_invalid_probability_inputs_fail_explicitly(matrix, counts):
    with pytest.raises(ValueError):
        validate_inputs(matrix, counts)


def test_zero_probability_is_infinite_only_for_observed_channels():
    matrix = np.array([[1.0], [0.0]])
    assert np.isinf(NegLogLik([1], matrix, [1, 1]))
    assert NegLogLik([1], matrix, [1, 0]) == 0
    with pytest.raises(ValueError, match="cannot explain"):
        fit_likelihood(matrix, [1, 1])


def test_tiny_positive_probabilities_keep_their_true_gradient():
    matrix = np.array([[1e-310, 2e-310], [1.0, 1.0]])
    weights = fit_likelihood(matrix, [100, 0])
    assert weights == pytest.approx([0, 1], abs=1e-6)


def test_background_can_supply_missing_known_signature_support():
    matrix = np.zeros((96, 1))
    matrix[0, 0] = 1
    counts = np.ones(96)
    with pytest.raises(ValueError, match="cannot explain"):
        decompose_mutational_profile_counts(counts, (matrix, ["A"]), enable_dummy=False)
    weights, _, _ = decompose_mutational_profile_counts(counts, (matrix, ["A"]), enable_dummy=True)
    background = np.array([x for _, x in get_dummy_signatures_lists()]).T
    assert np.column_stack([matrix, background]) @ weights == pytest.approx(counts / 96, abs=1e-6)


def test_unknown_method_fails_instead_of_selecting_frobenius():
    with pytest.raises(ValueError, match="Unknown fitting method"):
        decompose_mutational_profile_counts(
            COUNTS, (MATRIX, ["A", "B"]), func="typo", enable_dummy=False
        )


def test_python_and_cli_share_mle_default():
    parser = argparse.ArgumentParser()
    SignatureMenu(parser)
    assert parser.get_default("method") == DEFAULT_METHOD == "MLE"
    assert (
        inspect.signature(decompose_mutational_profile_counts).parameters["func"].default
        == DEFAULT_METHOD
    )
    default, _, _ = decompose_mutational_profile_counts(
        COUNTS, (MATRIX, ["A", "B"]), enable_dummy=False
    )
    explicit, _, _ = decompose_mutational_profile_counts(
        COUNTS, (MATRIX, ["A", "B"]), func="MLE", enable_dummy=False
    )
    assert np.array_equal(default, explicit)
