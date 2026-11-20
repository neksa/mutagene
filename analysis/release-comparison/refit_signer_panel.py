"""Refit the release default on the archived fixed-catalog signeR panel."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mutagene.io.matrix import channel_labels
from mutagene.signatures.identify import decompose_mutational_profile_counts


def refit(archive):
    root = Path(__file__).resolve().parents[2]
    manifest = json.loads((root / "analysis/local-artifacts.json").read_text())["files"]
    base = Path("analysis/signature-benchmark-2026")
    fresh = base / "background-results/fresh"
    metrics = base / "signer-results/metrics.tsv"
    for path in [fresh / name for name in ("inputs.npz", "reference.tsv", "samples.tsv")] + [
        metrics
    ]:
        if (
            hashlib.sha256((archive / path).read_bytes()).hexdigest()
            != manifest[str(path)]["sha256"]
        ):
            raise ValueError(f"Unexpected artifact checksum for {path}")

    with np.load(archive / fresh / "inputs.npz") as data:
        reference = data["reference"]
        names = data["names"].tolist()
        positions = [data["generating_names"].tolist().index(name) for name in names]
        counts = data["counts"]
        truth = data["truth"][:, positions]
    labelled = pd.read_csv(archive / fresh / "reference.tsv", sep="\t", index_col=0)
    assert labelled.index.tolist() == channel_labels() and labelled.columns.tolist() == names
    np.testing.assert_allclose(reference, labelled.to_numpy(), atol=1e-14)
    samples = pd.read_csv(archive / fresh / "samples.tsv", sep="\t")
    selected = samples.loc[
        samples.sample_index == 0, ["sample", "burden", "condition", "null_truth"]
    ]
    historical = pd.read_csv(
        archive / metrics, sep="\t", usecols=["sample", "method", "aggregate_null_tae"]
    ).pivot(index="sample", columns="method", values="aggregate_null_tae")
    assert len(selected) == 96 and set(historical.index) == set(selected["sample"])
    errors = []
    for row in selected.itertuples():
        weights, _, _ = decompose_mutational_profile_counts(counts[row.Index], (reference, names))
        known_error = abs(weights[: len(names)] - truth[row.Index]).sum()
        background_error = abs(weights[len(names) :].sum() - row.null_truth)
        errors.append((known_error + background_error) / 2)
    result = selected.assign(MutaGene=errors).set_index("sample").join(historical, validate="1:1")
    assert result.notna().all().all()
    high = result.burden > 1000
    np.testing.assert_allclose(
        result.loc[high, "MutaGene"], result.loc[high, "joint-BIC"], atol=1e-9
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="Root of the restored research files")
    parser.add_argument("output", type=Path, help="New TSV file for per-profile errors")
    args = parser.parse_args()
    refit(args.archive).to_csv(args.output, sep="\t", mode="x")
