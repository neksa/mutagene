"""Test subset selection against known mixtures and explicit background components."""

import itertools

import numpy as np
import pytest

from mutagene.signatures.sparse import fit_sparse, information_criterion


@pytest.mark.parametrize("criterion", ["aicc", "bic"])
def test_identifiable_mixture_retains_true_signatures(criterion):
    matrix = np.array([[0.8, 0.05, 0.05], [0.1, 0.85, 0.05], [0.1, 0.1, 0.9]])
    counts = 10000 * (matrix @ np.array([0.7, 0.3, 0]))
    info = {}
    weights = fit_sparse(matrix, counts, 3, criterion, info)
    assert np.allclose(weights, [0.7, 0.3, 0], atol=1e-5)
    assert info["selected_known_indices"] == [0, 1]
    assert info["optimizer_success"]
    values = [step["criterion"] for step in info["selection_trace"]]
    assert all(b < a for a, b in itertools.pairwise(values))


def test_background_can_explain_a_profile_without_known_signatures():
    matrix = np.array([[0.95, 0.05, 0.5], [0.05, 0.95, 0.5]])
    info = {}
    weights = fit_sparse(matrix, [500, 500], 2, diagnostics=info)
    assert np.allclose(weights, [0, 0, 1])
    assert info["selected_known_indices"] == []


def test_correlated_duplicates_do_not_both_survive():
    matrix = np.array([[0.9, 0.9, 0.1], [0.1, 0.1, 0.9]])
    weights = fit_sparse(matrix, [900, 100], 3)
    assert np.count_nonzero(weights > 1e-6) == 1
    assert np.allclose(matrix @ weights, [0.9, 0.1], atol=1e-5)


def test_single_signature_without_background_is_valid():
    assert np.allclose(fit_sparse(np.array([[0.8], [0.2]]), [80, 20], 1), [1])


def test_unexplained_observed_channels_fail_explicitly():
    with pytest.raises(ValueError, match="cannot explain"):
        fit_sparse(np.array([[1.0], [0.0]]), [80, 20], 1)


def test_aicc_respects_small_sample_limit():
    assert np.isinf(information_criterion(10, 5, 6, "aicc"))
    assert information_criterion(10, 1, 100, "bic") == pytest.approx(20 + np.log(100))


def test_failed_refits_do_not_silently_use_nnls(monkeypatch):
    from types import SimpleNamespace

    monkeypatch.setattr(
        "mutagene.signatures.sparse.minimize",
        lambda *a, **k: SimpleNamespace(
            success=False, nit=1, message="deliberate failure", x=np.array([0.5, 0.5])
        ),
    )
    with pytest.raises(RuntimeError, match="deliberate failure"):
        fit_sparse(np.array([[0.8, 0.2], [0.2, 0.8]]), [80, 20], 2)
