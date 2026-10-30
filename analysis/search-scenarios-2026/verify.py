"""Independently verify metric arithmetic, source identity, and transition endpoints."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]


def main():
    counts = {}
    for phase, repeats in [("screen", 6), ("transition", 3)]:
        folder = ROOT / phase
        data = dict(np.load(folder / "inputs.npz"))
        samples = pd.read_csv(folder / "samples.tsv", sep="\t")
        metrics = pd.read_csv(folder / "metrics.tsv", sep="\t")
        weights = np.load(folder / "estimates.npz")["weights"]
        assert len(metrics) == len(weights) == len(samples) * repeats
        assert (
            metrics["sample"].tolist() == np.repeat(samples["sample"].to_numpy(), repeats).tolist()
        )
        np.testing.assert_array_equal(data["counts"].sum(axis=1), samples.burden)
        np.testing.assert_array_equal(metrics.burden, np.repeat(samples.burden.to_numpy(), repeats))
        np.testing.assert_array_equal(
            metrics.null_truth, np.repeat(samples.null_truth.to_numpy(), repeats)
        )
        np.testing.assert_allclose(weights.sum(axis=1), 1, atol=1e-12, rtol=0)
        assert np.isfinite(weights).all() and (weights >= 0).all()
        k = len(data["names"])
        positions = [data["generating_names"].tolist().index(name) for name in data["names"]]
        truth = np.repeat(data["truth"], repeats, axis=0)
        known = weights[:, :k]
        null = weights[:, k:].sum(axis=1)
        expected = (
            np.abs(known - truth[:, positions]).sum(axis=1)
            + np.abs(null - metrics.null_truth.to_numpy())
        ) / 2
        np.testing.assert_allclose(expected, metrics.aggregate_null_tae, atol=1e-12, rtol=0)
        thresholded = known.copy()
        thresholded[known * metrics.burden.to_numpy()[:, None] < 10] = 0
        expanded = np.zeros_like(truth)
        expanded[:, positions] = thresholded
        np.testing.assert_allclose(
            np.abs(expanded - truth).sum(axis=1) / 2, metrics.paper_tae, atol=1e-12, rtol=0
        )
        spec = json.loads((folder / "specification.json").read_text())
        for name, digest in spec["sources"].items():
            assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == digest, name
        if phase == "transition":
            w = weights.reshape(len(samples), 3, -1)
            n = samples.burden.to_numpy()
            endpoints = (n <= 1000) | (n >= 2000)
            np.testing.assert_array_equal(w[endpoints, 0], w[endpoints, 1])
            middle = ~endpoints
            alpha = (n[middle] - 1000) / 1000
            np.testing.assert_allclose(
                w[middle, 1],
                alpha[:, None] * w[middle, 0] + (1 - alpha[:, None]) * w[middle, 2],
                atol=1e-12,
                rtol=0,
            )
        counts[phase] = len(metrics)
    assert (ROOT / "frozen_sparse.py").read_bytes() == (
        REPO / "mutagene/signatures/sparse.py"
    ).read_bytes()
    manifest = json.loads((ROOT / "artifact-manifest.json").read_text())
    for name, record in manifest.items():
        path = REPO / name
        assert path.stat().st_size == record["bytes"], name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"], name
    print(
        json.dumps(
            dict(
                metric_rows=counts,
                manifest_entries=len(manifest),
                production_default_unchanged=True,
                status="pass",
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
