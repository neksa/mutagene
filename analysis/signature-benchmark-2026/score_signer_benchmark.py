"""Score all comparators on the same signeR subset and verify its provenance."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from run_benchmark import score, sha256
from signer_adapter import read_fit

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "signer-results"
SOURCE = ROOT / "background-results/fresh"
specification = json.loads((OUT / "specification.json").read_text())
assert sha256(ROOT / "run_signer.R") == specification["adapter_sha256"]
assert sha256(ROOT / "prepare_signer_benchmark.py") == specification["preparation_sha256"]
for name, digest in specification["input_sha256"].items():
    assert sha256(SOURCE / name) == digest
inputs = dict(np.load(SOURCE / "inputs.npz"))
samples = pd.read_csv(OUT / "samples.tsv", sep="\t")
tasks = pd.read_csv(OUT / "tasks.tsv", sep="\t")
reference = pd.read_csv(OUT / "reference.tsv", sep="\t", index_col=0)
background = pd.read_csv(OUT / "background.tsv", sep="\t", index_col=0)
assert len(samples) == 96 and len(tasks) == 216
np.testing.assert_allclose(reference.to_numpy(), inputs["reference"], atol=1e-15, rtol=0)
catalogs = pd.read_csv(OUT / "catalogs.tsv", sep="\t", index_col=0)
assert catalogs.columns.tolist() == samples["sample"].tolist()
np.testing.assert_array_equal(catalogs.to_numpy().T, inputs["counts"][samples.input_index])
lookup = samples.set_index("sample")
rows = []
repeat_rows = []
knowns = []
nulls = []
for task in tasks.to_dict("records"):
    weights, native, info = read_fit(
        OUT / "fits" / f"{task['task_id']:04d}", task, reference.columns, background.columns
    )
    metadata = lookup.loc[task["sample"]].to_dict()
    i = int(metadata["input_index"])
    burden = inputs["counts"][i].sum()
    assert np.isclose(info["raw_total_ratio"], native.sum() / burden)
    known = weights[: reference.shape[1]]
    null = weights[reference.shape[1] :] if task["background"] else np.zeros(6)
    metrics = score(
        known,
        null,
        inputs["truth"][i],
        inputs["counts"][i],
        inputs["names"].tolist(),
        inputs["generating_names"].tolist(),
        inputs["reference"],
        metadata["null_truth"],
    )
    raw_known = native[: reference.shape[1]] / burden
    raw_null = native[reference.shape[1] :] / burden if task["background"] else np.zeros(6)
    raw_metrics = score(
        raw_known,
        raw_null,
        inputs["truth"][i],
        inputs["counts"][i],
        inputs["names"].tolist(),
        inputs["generating_names"].tolist(),
        inputs["reference"],
        metadata["null_truth"],
    )
    row = dict(
        metadata,
        sample=task["sample"],
        method="signeR+background" if task["background"] else "signeR",
        background=task["background"],
        chain=task["chain"],
        seed=task["seed"],
        seconds=info["seconds"],
        error=None,
        **metrics,
    )
    row["active_known"] = int(np.count_nonzero(known * burden >= 10))
    for key in ["raw_total_ratio", "active_split_rhat_max", "active_ess_bulk_min", "mcse_mean_max"]:
        row[key] = info[key]
    row["raw_aggregate_null_tae"] = raw_metrics["aggregate_null_tae"]
    row["raw_paper_tae"] = raw_metrics["paper_tae"]
    if task["chain"] == 1:
        rows.append(row)
        knowns.append(known)
        nulls.append(null)
    else:
        repeat_rows.append(row)
new = pd.DataFrame(rows)
new.to_csv(OUT / "signer-metrics.tsv", sep="\t", index=False)
pd.DataFrame(repeat_rows).to_csv(OUT / "repeat-metrics.tsv", sep="\t", index=False)
np.savez_compressed(OUT / "estimates.npz", known=np.array(knowns), null=np.array(nulls))
assert len(new) == 192 and not new.duplicated(["sample", "method"]).any()
baseline = pd.read_csv(SOURCE / "metrics.tsv", sep="\t", low_memory=False)
baseline = baseline[baseline["sample"].isin(samples["sample"])].copy()
assert len(baseline) == 576
assert set(baseline["sample"]) == set(new["sample"])
combined = pd.concat([baseline, new], ignore_index=True)
assert combined.groupby("method").size().eq(96).all()
combined.to_csv(OUT / "metrics.tsv", sep="\t", index=False)
metrics = [
    "aggregate_null_tae",
    "paper_tae",
    "active_known",
    "null_weight",
    "null_absolute_error",
    "cosine",
    "seconds",
]
combined.groupby(["condition", "burden", "method"])[metrics].mean().reset_index().to_csv(
    OUT / "summary.tsv", sep="\t", index=False
)
combined.groupby("method")[metrics].mean().reset_index().to_csv(
    OUT / "overall.tsv", sep="\t", index=False
)
combined.groupby(["cancer", "method"])[metrics].mean().reset_index().to_csv(
    OUT / "cancer-summary.tsv", sep="\t", index=False
)
# Resample cancer groups together; each group retains its twelve matched profiles.
rng = np.random.default_rng(202610281)
comparisons = []
for metric in ["aggregate_null_tae", "paper_tae"]:
    pivot = combined.pivot_table(
        index=["sample", "cancer", "condition", "burden"], columns="method", values=metric
    )
    for candidate in ["signeR", "signeR+background"]:
        for comparator in ["MLE", "fixed-BIC", "joint-BIC", "SigProfilerAssignment"]:
            diff = pivot[candidate] - pivot[comparator]
            cancer = diff.groupby(level="cancer").mean().to_numpy()
            draws = rng.choice(cancer, size=(2000, len(cancer)), replace=True).mean(axis=1)
            low, high = np.quantile(draws, [0.025, 0.975])
            comparisons.append(
                dict(
                    metric=metric,
                    candidate=candidate,
                    comparator=comparator,
                    difference=diff.mean(),
                    low=low,
                    high=high,
                )
            )
pd.DataFrame(comparisons).to_csv(OUT / "paired-differences.tsv", sep="\t", index=False)
checks = {
    "profiles": 96,
    "native_signeR_fits": 216,
    "primary_signeR_fits": 192,
    "independent_repeat_fits": 24,
    "matched_comparison_rows": 768,
    "errors": 0,
    "input_hashes": "pass",
    "signature_alignment": "pass",
    "relative_exposure_arithmetic": "pass",
    "raw_exposures_retained": True,
    "baseline_rows_unchanged": True,
    "caution": "Numerical completion does not certify MCMC convergence. Read chain-stability.tsv and per-fit diagnostics.",
}
(OUT / "verification.json").write_text(json.dumps(checks, indent=2) + "\n")
print(combined.groupby("method")[metrics].mean().to_string())
