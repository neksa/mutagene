"""Check benchmark accounting and the interpretation of omitted signature mass."""

import importlib.util
import sys
from pathlib import Path

import numpy as np


def load_runner():
    directory = Path(__file__).resolve().parents[2] / "analysis/signature-benchmark-2026"
    sys.path.insert(0, str(directory))
    try:
        spec = importlib.util.spec_from_file_location(
            "adoption_benchmark", directory / "run_benchmark.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.remove(str(directory))


def test_null_mass_does_not_identify_an_omitted_signature():
    runner = load_runner()
    matrix = np.full((96, 1), 1 / 96)
    counts = np.ones(96, dtype=int) * 100
    metrics = runner.score(
        np.array([0.8]),
        np.array([0.2, 0, 0, 0, 0, 0]),
        np.array([0.8, 0.2]),
        counts,
        ["known"],
        ["known", "omitted"],
        matrix,
        0.2,
    )
    assert metrics["aggregate_null_tae"] == 0
    assert np.isclose(metrics["paper_tae"], 0.1)
    assert metrics["known_weight"] == 0.8


def test_ten_mutation_cutoff_does_not_renormalize_remaining_weights():
    runner = load_runner()
    matrix = np.full((96, 2), 1 / 96)
    counts = np.zeros(96, dtype=int)
    counts[0] = 100
    metrics = runner.score(
        np.array([0.91, 0.09]),
        np.zeros(6),
        np.array([1.0, 0.0]),
        counts,
        ["a", "b"],
        ["a", "b"],
        matrix,
        0.0,
    )
    assert np.isclose(metrics["paper_tae"], 0.045)
    assert np.isclose(metrics["aggregate_null_tae"], 0.09)
    assert metrics["false_positive_weight"] == 0
