"""Verify repaired benchmark provenance, probability constraints, and metric arithmetic."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from run_benchmark import score, sha256

from mutagene.signatures import get_dummy_signatures_lists
from mutagene.signatures.likelihood import OPTIMALITY_TOLERANCE

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "repaired-results"
checks = {}
for phase, expected in [("paired", 20736), ("fresh", 22272)]:
    directory = OUT / phase
    specification = json.loads((directory / "specification.json").read_text())
    for filename, digest in specification["source_sha256"].items():
        assert sha256(ROOT.parents[1] / "mutagene/signatures" / filename) == digest
    assert sha256(ROOT / "run_repaired_benchmark.py") == specification["runner_sha256"]
    assert (
        sha256(ROOT / "baseline-20261026/mutagene/signatures/identify.py")
        == specification["baseline_source_sha256"]
    )
    with np.load(directory / "inputs.npz") as archive:
        inputs = {name: archive[name] for name in archive.files}
    with np.load(directory / "estimates.npz") as archive:
        estimates = {name: archive[name] for name in archive.files}
    metrics = pd.read_csv(directory / "metrics.tsv", sep="\t", low_memory=False)
    samples = pd.read_csv(directory / "samples.tsv", sep="\t")
    assert len(metrics) == expected and len(samples) == 768
    assert not metrics.duplicated(["sample", "method", "background"]).any()
    assert metrics.error.isna().all()
    assert np.array_equal(inputs["counts"].sum(axis=1), samples.burden)
    assert np.isfinite(estimates["known"]).all() and np.isfinite(estimates["null"]).all()
    assert (estimates["known"] >= 0).all() and (estimates["null"] >= 0).all()
    new = ~metrics.method.isin(["SigProfilerAssignment", "pre-repair-MLE"])
    sums = estimates["known"].sum(axis=1) + estimates["null"].sum(axis=1)
    assert np.allclose(sums[new], 1, atol=1e-12, rtol=0)
    likelihood = metrics.method.isin(["MLE", "MLE-scaled", "Compat"])
    sparse = metrics.method.str.startswith("sparse-")
    assert metrics.loc[likelihood | sparse, "fallback"].isna().all()
    assert (metrics.loc[likelihood, "optimality_gap"] <= OPTIMALITY_TOLERANCE).all()
    assert (metrics.loc[sparse, "maximum_refit_optimality_gap"] <= OPTIMALITY_TOLERANCE).all()
    indices = dict(zip(samples["sample"], samples.index))
    dummy = np.array([values for _, values in get_dummy_signatures_lists()]).T
    full = np.column_stack([inputs["reference"], dummy])
    # Check every repaired full-MLE certificate independently of the solver implementation.
    for index in metrics.index[likelihood]:
        row = metrics.iloc[index]
        i = indices[row["sample"]]
        weights = np.r_[estimates["known"][index], estimates["null"][index]]
        counts = inputs["counts"][i]
        observed = counts > 0
        design = full if row.background else full[:, :67]
        fitted = weights if row.background else weights[:67]
        predicted = design[observed] @ fitted
        p = counts[observed] / counts.sum()
        gap = max(0.0, (design[observed].T @ (p / predicted)).max() - 1)
        assert np.isclose(gap, row.optimality_gap, atol=1e-9, rtol=1e-10)
        ll = counts[observed] @ np.log(predicted)
        assert np.isclose(ll, row.count_log_likelihood, atol=1e-8, rtol=1e-12)
    for index in np.random.default_rng(20261027).choice(len(metrics), 200, replace=False):
        row = metrics.iloc[index]
        i = indices[row["sample"]]
        recomputed = score(
            estimates["known"][index],
            estimates["null"][index],
            inputs["truth"][i],
            inputs["counts"][i],
            inputs["names"].tolist(),
            inputs["generating_names"].tolist(),
            inputs["reference"],
            row.null_truth,
        )
        for key, value in recomputed.items():
            assert np.isclose(value, row[key], atol=1e-10, rtol=1e-10), (phase, index, key)
    if phase == "paired":
        with np.load(ROOT / "sparse-results/evaluate/inputs.npz") as original:
            for name in inputs:
                assert np.array_equal(inputs[name], original[name]), name
    checks[phase] = {
        "profiles": len(samples),
        "fits": len(metrics),
        "errors": 0,
        "certified_full_mle_fits": int(likelihood.sum()),
        "certified_sparse_fits": int(sparse.sum()),
        "metric_recomputations": 200,
        "maximum_mle_gap": float(metrics.loc[likelihood, "optimality_gap"].max()),
        "maximum_sparse_refit_gap": float(
            metrics.loc[sparse, "maximum_refit_optimality_gap"].max()
        ),
    }
    print(phase, checks[phase])
(OUT / "verification.json").write_text(json.dumps(checks, indent=2) + "\n")
