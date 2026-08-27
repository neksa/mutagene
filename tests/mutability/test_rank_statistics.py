"""Recurrence counting in rank.

The analysed sample adds a pseudocount to every observed mutation, but the
cohort size it is tested against was left unchanged. A mutation present in every
cohort sample then reached observed = N + 1, which predict_driver treats as
impossible, so the mutations with the strongest recurrence signal came out
labelled "Undefined".
"""

import pytest

from mutagene.mutability.mutability import predict_driver, rank


class TestPredictDriver:
    def test_a_mutation_in_every_sample_is_a_driver(self):
        _p, label, *_ = predict_driver(100, 100, 1e-7, 8.03e-05, 0.0034)
        assert label == "Driver"

    def test_observed_above_n_is_rejected_as_impossible(self):
        """This is the state the missing +1 produced."""
        pvalue, label, *_ = predict_driver(101, 100, 1e-7, 8.03e-05, 0.0034)
        assert (pvalue, label) == (1.0, "Undefined")

    def test_the_pseudocount_fits_once_the_cohort_counts_the_sample(self):
        _p, label, *_ = predict_driver(101, 101, 1e-7, 8.03e-05, 0.0034)
        assert label == "Driver"


class TestRankWithCohort:
    """rank() must not lose a mutation that every cohort sample carries."""

    def mutations(self):
        # GAA -> a glutamate codon; the substitution itself does not matter here
        return {("TP53", "R248W"): {"seq5": {"ACGGA": 1}, "transcript": "T1"}}

    def profile(self):
        return [10.0] * 96

    def test_ubiquitous_mutation_is_not_undefined(self, tmp_path):
        out = tmp_path / "out.tsv"
        cohort_size = 50
        cohort = {"TP53": {"R248W": cohort_size}}  # present in every sample

        with open(out, "w") as fh:
            rank(self.mutations(), fh, self.profile(), cohort, cohort_size, 8.03e-05, 0.0034)

        rows = [line for line in out.read_text().splitlines() if not line.startswith("#")]
        assert len(rows) >= 2, "the mutation was dropped entirely"
        assert "Undefined" not in rows[1], rows[1]

    def test_a_rare_mutation_is_unaffected(self, tmp_path):
        out = tmp_path / "out.tsv"
        cohort = {"TP53": {"R248W": 1}}

        with open(out, "w") as fh:
            rank(self.mutations(), fh, self.profile(), cohort, 50, 8.03e-05, 0.0034)

        rows = [line for line in out.read_text().splitlines() if not line.startswith("#")]
        assert len(rows) >= 2


def test_cohort_size_is_untouched_without_a_cohort(tmp_path):
    """Ranking against the input sample alone adds no pseudocount, so no +1."""
    out = tmp_path / "out.tsv"
    mutations = {("TP53", "R248W"): {"seq5": {"ACGGA": 1}, "transcript": "T1"}}

    with open(out, "w") as fh:
        rank(mutations, fh, [10.0] * 96, None, 1, 8.03e-05, 0.0034)

    assert out.exists()
