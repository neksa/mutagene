"""Independently recalculate all validation errors and verify the final artifact hashes."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]


def main():
    counts = {}
    for phase in ["validation", "sensitivity"]:
        folder = ROOT / phase
        inputs = dict(np.load(folder / "inputs.npz"))
        samples = pd.read_csv(folder / "samples.tsv", sep="\t")
        metrics = pd.read_csv(folder / "metrics.tsv", sep="\t")
        weights = np.load(folder / "estimates.npz")["weights"]
        assert len(metrics) == len(samples) * 2 == len(weights)
        assert metrics["sample"].tolist() == np.repeat(samples["sample"].to_numpy(), 2).tolist()
        assert metrics.method.tolist() == ["joint-BIC", "adaptive-BIC"] * len(samples)
        np.testing.assert_array_equal(inputs["counts"].sum(axis=1), samples.burden)
        np.testing.assert_array_equal(metrics.burden, np.repeat(samples.burden.to_numpy(), 2))
        np.testing.assert_array_equal(
            metrics.null_truth, np.repeat(samples.null_truth.to_numpy(), 2)
        )
        known_count = len(inputs["names"])
        positions = [inputs["generating_names"].tolist().index(name) for name in inputs["names"]]
        truth = np.repeat(inputs["truth"], 2, axis=0)
        known = weights[:, :known_count]
        null = weights[:, known_count:].sum(axis=1)
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
        np.testing.assert_allclose(weights.sum(axis=1), 1, atol=1e-12, rtol=0)
        assert (weights >= 0).all() and np.isfinite(weights).all()
        paired = weights.reshape(len(samples), 2, -1)
        high = samples.burden.to_numpy() > 1000
        np.testing.assert_array_equal(paired[high, 0], paired[high, 1])
        spec = json.loads((folder / "specification.json").read_text())
        for name, digest in spec["source_sha256"].items():
            assert (
                hashlib.sha256((folder / "source" / name).read_bytes()).hexdigest() == digest
            ), name
        assert (
            hashlib.sha256((ROOT / "selection.json").read_bytes()).hexdigest()
            == spec["selection_sha256"]
        )
        for name, digest in spec["package_source_sha256"].items():
            source = (
                ROOT / "baseline-source/mutagene/signatures" / name
                if name == "sparse.py"
                else REPO / "mutagene/signatures" / name
            )
            assert hashlib.sha256(source.read_bytes()).hexdigest() == digest, name
        counts[phase] = len(metrics)
    manifest = json.loads((ROOT / "artifact-manifest.json").read_text())
    for name, record in manifest.items():
        path = REPO / name
        assert path.stat().st_size == record["bytes"], name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"], name
    print(
        json.dumps(
            dict(metric_rows_verified=counts, manifest_entries=len(manifest), status="pass"),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
