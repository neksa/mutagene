"""Verify likelihood audit provenance, feasibility, and independent metric arithmetic."""

import json

import numpy as np
import pandas as pd
from audit_likelihood import OUT, ROOT
from run_benchmark import sha256

from mutagene.signatures import get_dummy_signatures_lists

provenance = json.loads((OUT / "provenance.json").read_text())
assert sha256(ROOT / "audit_likelihood.py") == provenance["audit_script_sha256"]
for name, digest in provenance["source_hashes"].items():
    assert sha256(ROOT.parents[1] / "mutagene/signatures" / name) == digest
source = ROOT / "sparse-results/evaluate"
for name, digest in provenance["input_hashes"].items():
    assert sha256(source / name) == digest
with np.load(source / "inputs.npz") as archive:
    inputs = {name: archive[name] for name in archive.files}
weights = np.load(OUT / "stable-estimates.npz")["weights"]
metrics = pd.read_csv(OUT / "metrics.tsv", sep="\t", low_memory=False)
stable = metrics[metrics.method == "stable-full-MLE"].reset_index(drop=True)
samples = pd.read_csv(source / "samples.tsv", sep="\t")
assert weights.shape == (1536, 73)
assert stable.success.all()
assert np.isfinite(weights).all() and (weights >= 0).all()
assert np.allclose(weights.sum(axis=1), 1, atol=1e-12)
assert stable.gap_bound.max() < 2e-5
assert (metrics.loc[metrics.method == "MLE", "nll_excess_over_stable"] > -1e-10).all()
dummy = np.array([values for _, values in get_dummy_signatures_lists()]).T
matrix = np.column_stack([inputs["reference"], dummy])
positions = [inputs["generating_names"].tolist().index(name) for name in inputs["names"]]
for j, row in stable.iterrows():
    i = j // 2
    assert row["sample"] == samples.iloc[i]["sample"]
    assert bool(row.background) == bool(j % 2)
    counts = inputs["counts"][i]
    positive = counts > 0
    q = matrix @ weights[j]
    nll = -np.dot(counts[positive] / counts.sum(), np.log(q[positive]))
    assert np.isclose(nll, row.nll_per_mutation, atol=1e-12)
    target = inputs["truth"][i, positions]
    error = (
        np.abs(weights[j, :67] - target).sum()
        + abs(weights[j, 67:].sum() - samples.iloc[i].null_truth)
    ) / 2
    assert np.isclose(error, row.aggregate_null_tae, atol=1e-12)
    active_matrix = matrix if row.background else matrix[:, :67]
    gap = (active_matrix[positive].T @ (counts[positive] / counts.sum() / q[positive])).max() - 1
    assert np.isclose(gap, row.gap_bound, atol=1e-10)
failures = pd.read_csv(OUT / "discarded-candidates.tsv", sep="\t")
assert len(failures) == 332 and not failures.optimizer_success.any()
assert (failures.fallback_nll_increase > 1e-10).sum() == 327
validation = {
    "source_and_input_hashes": "pass",
    "stable_simplex_constraints": "pass",
    "independent_nll_attribution_and_optimality_checks": 1536,
    "reproduced_optimizer_failures": 332,
    "normalized_candidates_better_than_fallback": 327,
}
(OUT / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
print(json.dumps(validation, indent=2))
