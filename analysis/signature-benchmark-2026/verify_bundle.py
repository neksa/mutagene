"""Verify frozen inputs, metric arithmetic, and comparator channel alignment."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from run_benchmark import channel_order, score, sha256

root = Path(__file__).resolve().parent
results = root / "results"
provenance = json.loads((results / "provenance.json").read_text())
assert sha256(root / "frozen_identify.py") == provenance["fitter_sha256"]
for filename, expected in provenance["output_hashes"].items():
    assert sha256(results / filename) == expected, filename
inputs = np.load(results / "inputs.npz")
estimates = np.load(results / "estimates.npz")
metrics = pd.read_csv(results / "metrics.tsv", sep="\t")
metadata = pd.read_csv(results / "samples.tsv", sep="\t")
reference = pd.read_csv(results / "reference.tsv", sep="\t", index_col=0)
assert len(metadata) == 4800
assert len(metrics) == 24000
assert len(metrics[["sample", "method"]].drop_duplicates()) == len(metrics)
assert np.array_equal(inputs["counts"].sum(axis=1), metadata.burden)
assert np.allclose(reference.sum(), 1)
assert reference.index.tolist() == channel_order()
assert np.isfinite(estimates["known"]).all()
assert np.isfinite(estimates["null"]).all()
positions = {name: i for i, name in enumerate(metadata["sample"])}
for index in np.random.default_rng(42).choice(len(metrics), 100, replace=False):
    row = metrics.iloc[index]
    sample_index = positions[row["sample"]]
    recomputed = score(
        estimates["known"][index],
        estimates["null"][index],
        inputs["truth"][sample_index],
        inputs["counts"][sample_index],
        inputs["signature_names"].tolist(),
        inputs["generating_names"].tolist(),
        inputs["reference"],
        row.null_truth,
    )
    for key, value in recomputed.items():
        assert np.isclose(value, row[key]), (index, key)
spa_reference = pd.read_csv(results / "spa_reference.tsv", sep="\t", index_col=0)
assert np.allclose(spa_reference.loc[reference.index], reference)
print(
    "Verified 4,800 inputs, 24,000 fit rows, reference alignment, hashes, and 100 metric recomputations"
)
