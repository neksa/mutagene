"""Check background benchmark arithmetic and constrained refit certificates."""

import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
from run_benchmark import score, sha256

from mutagene.signatures import get_dummy_signatures_lists
from mutagene.signatures.likelihood import OPTIMALITY_TOLERANCE

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "background-results"
checks = {}
for phase in ("paired", "fresh"):
    directory = OUT / phase
    specification = json.loads((directory / "specification.json").read_text())
    for filename, digest in specification["source_sha256"].items():
        assert sha256(OUT / "source" / filename) == digest
    assert sha256(OUT / "source/run_background_benchmark.py") == specification["runner_sha256"]
    inputs = dict(np.load(directory / "inputs.npz"))
    estimates = dict(np.load(directory / "estimates.npz"))
    metrics = pd.read_csv(directory / "metrics.tsv", sep="\t", low_memory=False)
    samples = pd.read_csv(directory / "samples.tsv", sep="\t")
    assert len(metrics) == 4608 and len(samples) == 768
    assert metrics.error.isna().all()
    assert not metrics.duplicated(["sample", "method"]).any()
    weights = np.column_stack([estimates["known"], estimates["null"]])
    assert np.isfinite(weights).all() and (weights >= 0).all()
    ours = metrics.method != "SigProfilerAssignment"
    assert np.allclose(weights[ours].sum(axis=1), 1, atol=1e-12, rtol=0)
    assert metrics.loc[ours, "fallback"].isna().all()
    sparse = ours & (metrics.method != "MLE")
    assert (metrics.loc[sparse, "maximum_refit_optimality_gap"] <= OPTIMALITY_TOLERANCE).all()
    indices = dict(zip(samples["sample"], samples.index))
    dummy = np.array([v for _, v in get_dummy_signatures_lists()]).T
    full = np.column_stack([inputs["reference"], dummy])
    full /= full.sum(axis=0)
    known_count = inputs["reference"].shape[1]
    for index in metrics.index[ours]:
        row = metrics.iloc[index]
        i = indices[row["sample"]]
        counts = inputs["counts"][i]
        if row["method"] == "MLE":
            columns = list(range(full.shape[1]))
        else:
            columns = ast.literal_eval(row.selected_known_indices) + [
                known_count + j for j in ast.literal_eval(row.selected_background_indices)
            ]
            assert row.model_dimensions == len(columns) - 1
        observed = counts > 0
        design = full[observed][:, columns]
        fitted = weights[index, columns]
        predicted = design @ fitted
        assert (predicted > 0).all()
        gap = max(0.0, (design.T @ (counts[observed] / counts.sum() / predicted)).max() - 1)
        assert gap <= OPTIMALITY_TOLERANCE + 1e-10
        if row["method"] != "MLE":
            bic = -2 * np.dot(counts[observed], np.log(predicted)) + (len(columns) - 1) * np.log(
                counts.sum()
            )
            assert np.isclose(bic, row.criterion_value, rtol=1e-11, atol=1e-7)
    for index in np.random.default_rng(202610274).choice(len(metrics), 200, replace=False):
        row = metrics.iloc[index]
        i = indices[row["sample"]]
        expected = score(
            estimates["known"][index],
            estimates["null"][index],
            inputs["truth"][i],
            inputs["counts"][i],
            inputs["names"].tolist(),
            inputs["generating_names"].tolist(),
            inputs["reference"],
            row.null_truth,
        )
        for key, value in expected.items():
            assert np.isclose(row[key], value, atol=1e-12)
    if phase == "paired":
        original = dict(np.load(ROOT / "repaired-results/fresh/inputs.npz"))
        for key in original:
            np.testing.assert_array_equal(inputs[key], original[key])
        prior = pd.read_csv(ROOT / "repaired-results/fresh/metrics.tsv", sep="\t", low_memory=False)
        for label, method, background in [
            ("MLE", "MLE", True),
            ("fixed-BIC", "sparse-bic", True),
            ("no-background-BIC", "sparse-bic", False),
            ("SigProfilerAssignment", "SigProfilerAssignment", False),
        ]:
            old = prior[(prior.method == method) & (prior.background == background)].set_index(
                "sample"
            )
            new = metrics[metrics.method == label].set_index("sample")
            np.testing.assert_allclose(
                new.paper_tae, old.loc[new.index].paper_tae, atol=1e-12, rtol=0
            )
    checks[phase] = {
        "fits": len(metrics),
        "errors": 0,
        "fallbacks": 0,
        "independent_certificates": int(ours.sum()),
        "metric_recomputations": 200,
        "maximum_refit_gap": float(metrics.loc[sparse, "maximum_refit_optimality_gap"].max()),
    }
(OUT / "verification.json").write_text(json.dumps(checks, indent=2) + "\n")
print(json.dumps(checks, indent=2))
