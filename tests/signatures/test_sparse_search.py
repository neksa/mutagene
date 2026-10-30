"""Protect low-burden search gains and preserve the high-burden fitting path."""

import json
from pathlib import Path

import numpy as np
import pytest

from mutagene.io.profile import read_signatures
from mutagene.signatures import get_dummy_signatures_lists
from mutagene.signatures.sparse import fit_sparse


def test_single_start_escapes_full_start_local_optimum():
    fixture = json.loads((Path(__file__).parent / "data/sparse-search-start.json").read_text())
    reference, names = read_signatures("COSMICv3")
    dummy = np.array([values for _, values in get_dummy_signatures_lists()]).T
    matrix = np.column_stack([reference, dummy])
    old, new = {}, {}
    fit_sparse(matrix, fixture["counts"], len(names), diagnostics=old, search_start="full")
    weights = fit_sparse(matrix, fixture["counts"], len(names), diagnostics=new)
    assert new["criterion_value"] < old["criterion_value"] - 10
    assert new["model_dimensions"] < old["model_dimensions"]
    assert new["search_start"] == "single"
    assert new["maximum_refit_optimality_gap"] <= 2e-5
    assert weights.sum() == pytest.approx(1)
    assert np.isfinite(weights).all() and (weights >= 0).all()


@pytest.mark.parametrize("burden", [1001, 2000, 50000])
def test_high_burden_weights_match_full_start_exactly(burden):
    matrix = np.array([[0.8, 0.1, 0.4], [0.15, 0.7, 0.4], [0.05, 0.2, 0.2]])
    counts = np.array([0.6, 0.3, 0.1]) * burden
    expected = fit_sparse(matrix, counts, 2, search_start="full")
    info = {}
    actual = fit_sparse(matrix, counts, 2, diagnostics=info)
    np.testing.assert_array_equal(actual, expected)
    assert info["search_start"] == "full"


def test_threshold_includes_one_thousand_mutations():
    info = {}
    fit_sparse(np.array([[0.8, 0.1], [0.2, 0.9]]), [800, 200], 2, diagnostics=info)
    assert info["search_start"] == "single"


@pytest.mark.parametrize("policy", ["fixed", "block"])
def test_other_background_policies_keep_full_start(policy):
    matrix = np.array([[0.8, 0.1, 0.5], [0.2, 0.9, 0.5]])
    expected = fit_sparse(matrix, [65, 35], 2, background_selection=policy, search_start="full")
    info = {}
    actual = fit_sparse(matrix, [65, 35], 2, background_selection=policy, diagnostics=info)
    np.testing.assert_array_equal(actual, expected)
    assert info["search_start"] == "full"


def test_aicc_keeps_full_start():
    matrix = np.array([[0.8, 0.1], [0.2, 0.9]])
    expected = fit_sparse(matrix, [65, 35], 2, criterion="aicc", search_start="full")
    info = {}
    actual = fit_sparse(matrix, [65, 35], 2, criterion="aicc", diagnostics=info)
    np.testing.assert_array_equal(actual, expected)
    assert info["search_start"] == "full"


def test_disjoint_support_uses_feasible_full_start():
    info = {}
    weights = fit_sparse(np.eye(3), [60, 30, 10], 2, diagnostics=info)
    np.testing.assert_allclose(weights, [0.6, 0.3, 0.1], atol=1e-7)
    assert info["search_start"] == "full"


def test_single_component_reference_finishes():
    info = {}
    weights = fit_sparse(np.array([[0.4], [0.6]]), [4, 6], 1, diagnostics=info)
    np.testing.assert_array_equal(weights, [1])
    assert info["model_dimensions"] == 0


def test_invalid_search_start_fails():
    with pytest.raises(ValueError, match="Search start"):
        fit_sparse(np.eye(2), [20, 80], 2, search_start="invalid")
