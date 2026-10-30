"""Check background selection, model dimensions, and complete BIC comparisons."""

import numpy as np
import pytest

from mutagene.signatures.sparse import fit_sparse, information_criterion


@pytest.mark.parametrize("policy", ["joint", "block"])
def test_clean_profile_removes_background(policy):
    matrix = np.array([[0.8, 0.1, 0.5], [0.2, 0.9, 0.5]])
    info = {}
    weights = fit_sparse(matrix, [8000, 2000], 2, diagnostics=info, background_selection=policy)
    np.testing.assert_allclose(weights, [1, 0, 0], atol=1e-6)
    assert info["selected_background_indices"] == []
    assert info["model_dimensions"] == 0


@pytest.mark.parametrize("policy", ["fixed", "joint", "block"])
def test_background_only_support_and_dimensions(policy):
    matrix = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    info = {}
    weights = fit_sparse(matrix, [0, 0, 1000], 1, diagnostics=info, background_selection=policy)
    np.testing.assert_allclose(weights, [0, 0, 1], atol=1e-12)
    expected = 0 if policy == "joint" else 1
    assert info["model_dimensions"] == expected
    assert info["criterion_value"] == pytest.approx(expected * np.log(1000))
    assert info["maximum_refit_optimality_gap"] <= 2e-5


def test_joint_retains_one_background_with_known_signature():
    matrix = np.eye(4)
    info = {}
    weights = fit_sparse(
        matrix, [8000, 0, 2000, 0], 2, diagnostics=info, background_selection="joint"
    )
    np.testing.assert_allclose(weights, [0.8, 0, 0.2, 0])
    assert info["selected_known_indices"] == [0]
    assert info["selected_background_indices"] == [0]
    assert info["model_dimensions"] == 1
    nll = -(8000 * np.log(0.8) + 2000 * np.log(0.2))
    assert info["criterion_value"] == pytest.approx(information_criterion(nll, 1, 10000, "bic"))


def test_block_matches_independent_complete_model_comparison():
    matrix = np.array([[0.8, 0.1, 0.5, 0.4], [0.2, 0.9, 0.5, 0.6]])
    counts = [650, 350]
    full, reduced, block = {}, {}, {}
    w_full = fit_sparse(matrix, counts, 2, diagnostics=full, background_selection="fixed")
    w_reduced = fit_sparse(matrix[:, :2], counts, 2, diagnostics=reduced)
    weights = fit_sparse(matrix, counts, 2, diagnostics=block, background_selection="block")
    best = min(full["criterion_value"], reduced["criterion_value"])
    assert block["criterion_value"] == best
    expected = w_full if best == full["criterion_value"] else np.r_[w_reduced, 0, 0]
    np.testing.assert_allclose(weights, expected)
    assert block["subset_refits"] == full["subset_refits"] + reduced["subset_refits"]


def test_invalid_policy_fails():
    with pytest.raises(ValueError, match="Background selection"):
        fit_sparse(np.eye(2), [10, 20], 1, background_selection="invalid")


@pytest.mark.parametrize("policy", ["joint", "block"])
def test_no_background_matches_fixed(policy):
    matrix = np.array([[0.8, 0.1], [0.2, 0.9]])
    expected = fit_sparse(matrix, [800, 200], 2)
    actual = fit_sparse(matrix, [800, 200], 2, background_selection=policy)
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("policy", ["joint", "block"])
def test_selection_is_reproducible(policy):
    matrix = np.array([[0.75, 0.10, 0.50, 0.25], [0.25, 0.90, 0.50, 0.75]])
    first, second = {}, {}
    a = fit_sparse(matrix, [650, 350], 2, diagnostics=first, background_selection=policy)
    b = fit_sparse(matrix, [650, 350], 2, diagnostics=second, background_selection=policy)
    np.testing.assert_array_equal(a, b)
    assert first == second


@pytest.mark.parametrize("policy", ["joint", "block"])
def test_candidate_failure_propagates(policy, monkeypatch):
    def fail(*args, **kwargs):
        raise RuntimeError("Deliberate refit failure")

    monkeypatch.setattr("mutagene.signatures.sparse.fit_likelihood", fail)
    with pytest.raises(RuntimeError, match="Deliberate refit failure"):
        fit_sparse(np.eye(2), [100, 50], 1, background_selection=policy)


def test_bic_selection_responds_to_mutation_burden():
    matrix = np.array([[0.8, 0.2], [0.2, 0.8]])
    low, high = {}, {}
    fit_sparse(matrix, [7, 3], 2, diagnostics=low)
    fit_sparse(matrix, [700, 300], 2, diagnostics=high)
    assert low["selected_known_indices"] == [0]
    assert high["selected_known_indices"] == [0, 1]
