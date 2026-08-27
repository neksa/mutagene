"""Signature decomposition: likelihood and channel ordering.

This module carries the arithmetic that produced wrong results before -- a
profile paired with the signature matrix by sorted key rather than by row order
gives a plausible-looking decomposition of the wrong data. Nothing tested it.
"""

import numpy as np
import pytest

from mutagene.io.profile import get_profile_attributes_dict, read_signatures
from mutagene.signatures.identify import NegLogLik, decompose_mutational_profile_counts
from mutagene.webapp.analysis import profile_channel_order, profile_dict_to_array


class TestChannelOrder:
    """The invariant documented in profile_doc.rst, now enforced."""

    def test_ninety_six_channels(self):
        assert len(profile_channel_order()) == 96
        assert len(set(profile_channel_order())) == 96

    def test_order_is_not_alphabetical(self):
        """Sorting the labels gives a different permutation of the same channels."""
        order = profile_channel_order()
        assert order != sorted(order)

    def test_labels_match_the_signature_matrix_rows(self):
        attributes = get_profile_attributes_dict()
        for label, attr in zip(profile_channel_order(), attributes):
            five, three = attr["context"]
            ref, alt = attr["mutation"]
            assert label == f"{five}[{ref}>{alt}]{three}"

    def test_dict_to_array_follows_that_order(self):
        order = profile_channel_order()
        counts = {label: i for i, label in enumerate(order)}
        assert np.array_equal(profile_dict_to_array(counts), np.arange(96, dtype=float))

    def test_a_signature_matrix_has_one_row_per_channel(self):
        W, names = read_signatures("COSMICv3", only=None)
        assert W.shape[0] == len(profile_channel_order())
        assert W.shape[1] == len(names)


class TestNegLogLik:
    def test_a_possible_mixture_scores_better_than_an_impossible_one(self):
        # channel 0 is observed, but signature 1 gives it zero probability
        A = np.array([[0.5, 0.0], [0.5, 1.0]])
        b = np.array([10.0, 10.0])

        possible = NegLogLik(np.array([1.0, 0.0]), A, b)
        impossible = NegLogLik(np.array([0.0, 1.0]), A, b)

        assert impossible > possible, "a mixture ruling out an observed channel scored no worse"

    def test_an_unobserved_impossible_channel_is_not_penalised(self):
        """Zero probability for something never seen is not evidence against."""
        A = np.array([[0.5, 0.0], [0.5, 1.0]])
        b = np.array([0.0, 10.0])

        assert NegLogLik(np.array([0.0, 1.0]), A, b) == pytest.approx(0.0, abs=1e-9)


class TestDecomposition:
    def test_a_pure_signature_is_recovered(self):
        W, names = read_signatures("COSMICv2", only=None)
        target = W[:, 0] * 5000  # 5000 mutations, all from signature 1

        _a, _b, results = decompose_mutational_profile_counts(
            target, (W, names), func="MLE", others_threshold=0.0, enable_dummy=False
        )

        best = max(results, key=lambda r: r.get("score", 0))
        assert best["name"] == names[0]
        assert best["score"] > 0.8, f"recovered only {best['score']:.2f} of a pure signature"

    def test_an_even_mixture_of_two_signatures_is_recovered(self):
        W, names = read_signatures("COSMICv2", only=None)
        target = (W[:, 0] + W[:, 4]) * 2500

        _a, _b, results = decompose_mutational_profile_counts(
            target, (W, names), func="MLE", others_threshold=0.0, enable_dummy=False
        )

        # The tail of the list holds diagnostics, not signatures.
        exposures = {r["name"]: r["score"] for r in results if r["name"] in set(names)}

        assert sum(exposures.values()) == pytest.approx(1.0, abs=0.05)
        assert exposures[names[0]] == pytest.approx(0.5, abs=0.1)
        assert exposures[names[4]] == pytest.approx(0.5, abs=0.1)

    def test_global_optimization_is_refused_rather_than_crashing(self):
        W, names = read_signatures("COSMICv2", only=None)

        with pytest.raises(NotImplementedError, match="global_optimization"):
            decompose_mutational_profile_counts(
                W[:, 0] * 100, (W, names), func="MLE", global_optimization=True
            )
