"""Audit sparse evaluation provenance, pairings, estimates, and metric arithmetic."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from run_benchmark import score, sha256

root = Path(__file__).resolve().parent
bundle = root / "sparse-results"
selection = json.loads((bundle / "selection.json").read_text())
assert sha256(bundle / "tune/metrics.tsv") == selection["tuning_metrics_sha256"]
assert (
    sha256(root.parents[1] / "mutagene/signatures/sparse.py") == selection["sparse_source_sha256"]
)
for phase, profiles, fits in [("tune", 192, 768), ("evaluate", 768, 20736)]:
    directory = bundle / phase
    spec = json.loads((directory / "specification.json").read_text())
    metrics = pd.read_csv(directory / "metrics.tsv", sep="\t", low_memory=False)
    metadata = pd.read_csv(directory / "samples.tsv", sep="\t")
    inputs = np.load(directory / "inputs.npz")
    truth, counts = inputs["truth"], inputs["counts"]
    estimates = np.load(directory / "estimates.npz")
    known, null = estimates["known"], estimates["null"]
    assert len(metrics) == fits and len(metadata) == profiles
    assert not metrics.duplicated(["sample", "method", "background"]).any()
    assert np.array_equal(counts.sum(axis=1), metadata.burden)
    assert metrics.error.isna().all(), "Unresolved errors need explicit reporting"
    assert np.isfinite(known).all() and np.isfinite(null).all()
    for filename, digest in spec["source_sha256"].items():
        assert sha256(root.parents[1] / "mutagene/signatures" / filename) == digest
    if phase == "evaluate":
        assert spec["frozen_selection"] == selection
    lookup = dict(zip(metadata["sample"], range(len(metadata))))
    for index in np.random.default_rng(42).choice(len(metrics), 100, replace=False):
        row = metrics.iloc[index]
        i = lookup[row["sample"]]
        recomputed = score(
            known[index],
            null[index],
            truth[i],
            counts[i],
            inputs["names"].tolist(),
            inputs["generating_names"].tolist(),
            inputs["reference"],
            row.null_truth,
        )
        for key, value in recomputed.items():
            assert np.isclose(value, row[key]), (phase, index, key)
        assert int(np.count_nonzero(known[index] * row.burden >= 10)) == row.active_known
    sparse = metrics.method.str.startswith("sparse-").to_numpy()
    assert np.allclose(known[sparse].sum(axis=1) + null[sparse].sum(axis=1), 1)
    assert (known[sparse] >= 0).all() and (null[sparse] >= 0).all()
    assert metrics.loc[sparse, "fallback"].isna().all()
print(
    "Verified tuning and evaluation provenance, all pairings, sparse simplex constraints, and 200 metric recomputations"
)
