import numpy as np
import pandas as pd
import pytest

from mutagene.signatures.discovery import compare_candidate


def test_disjoint_candidate_recovers_missing_channel(tmp_path):
    counts = pd.DataFrame(
        [[70, 0, 50, 20, 0, 85, 70, 0, 1], [30, 100, 50, 80, 100, 15, 30, 100, 99]],
        index=["A", "B"],
    )
    ref = pd.DataFrame({"known": [1, 0]}, index=counts.index)
    counts.to_csv(tmp_path / "counts.tsv", sep="\t")
    ref.iloc[::-1].to_csv(tmp_path / "reference.tsv", sep="\t")
    output = tmp_path / "results"
    summary = compare_candidate(tmp_path / "counts.tsv", tmp_path / "reference.tsv", output).iloc[0]
    assert summary.mixture_cosine == pytest.approx(0)
    assert summary.mean_gain > 0
    assert (summary.training_tumors, summary.test_tumors) == (6, 3)
    profile = pd.read_csv(output / "candidate-profiles.tsv", sep="\t", index_col=0)
    assert profile.loc["B", "Full catalog"] == pytest.approx(1)
    assert np.isnan(profile.loc["A", "Full catalog mixture"])


def test_exact_catalog_explains_training_profiles(tmp_path):
    counts = pd.DataFrame(np.ones((2, 6)))
    counts.to_csv(tmp_path / "counts.tsv", sep="\t")
    pd.DataFrame({"known": [1, 1]}).to_csv(tmp_path / "reference.tsv", sep="\t")
    with pytest.raises(ValueError, match="explain the training profiles"):
        compare_candidate(tmp_path / "counts.tsv", tmp_path / "reference.tsv", tmp_path)
