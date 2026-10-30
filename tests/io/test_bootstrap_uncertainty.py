"""Bootstrap uncertainty must not collapse when Monte Carlo replication increases."""

import numpy as np

from mutagene.io.decomposition import _get_bootstrap_stats_t


def test_repeated_bootstrap_draws_preserve_sampling_uncertainty():
    sample = [{"name": "SBS1", "score": 0.5, "mutations": 50}]
    draws = [
        [{"name": "SBS1", "score": x, "mutations": x * 100}] for x in [0.3, 0.4, 0.5, 0.6, 0.7]
    ]
    small = _get_bootstrap_stats_t(sample, draws, 100, 95).iloc[0]
    large = _get_bootstrap_stats_t(sample, draws * 100, 100, 95).iloc[0]
    width_small = small.exposure_high - small.exposure_low
    width_large = large.exposure_high - large.exposure_low
    assert 0.85 < width_large / width_small < 1.05
    assert np.isclose(large.exposure, 0.5)


def test_bias_correction_preserves_physical_bounds():
    for estimate, draw in [(0.0, 0.9), (1.0, 0.1)]:
        sample = [{"name": "SBS1", "score": estimate, "mutations": estimate * 100}]
        draws = [[{"name": "SBS1", "score": draw, "mutations": draw * 100}]] * 10
        result = _get_bootstrap_stats_t(sample, draws, 100, 95).iloc[0]
        assert 0 <= result.exposure_low <= result.exposure_high <= 1
        assert 0 <= result.mutations_low <= result.mutations_high <= 100
