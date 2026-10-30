"""Verify current defaults against the explicitly configured benchmark candidate."""

import inspect
import json

import numpy as np
import pandas as pd
from run_benchmark import ROOT as REPOSITORY_ROOT

from mutagene.signatures.identify import decompose_mutational_profile_counts
from mutagene.signatures.likelihood import DEFAULT_BACKGROUND_SELECTION, DEFAULT_METHOD
from mutagene.signatures.sparse import fit_sparse

ROOT = REPOSITORY_ROOT / "analysis/signature-benchmark-2026"
OUT = ROOT / "background-results"
decision = json.loads((OUT / "default-decision.json").read_text())
assert DEFAULT_METHOD == decision["method"]
assert DEFAULT_BACKGROUND_SELECTION == decision["background_selection"]
assert (
    inspect.signature(decompose_mutational_profile_counts).parameters["func"].default
    == DEFAULT_METHOD
)
assert (
    inspect.signature(decompose_mutational_profile_counts)
    .parameters["background_selection"]
    .default
    == DEFAULT_BACKGROUND_SELECTION
)
assert (
    inspect.signature(fit_sparse).parameters["background_selection"].default
    == DEFAULT_BACKGROUND_SELECTION
)
inputs = dict(np.load(OUT / "fresh/inputs.npz"))
estimates = dict(np.load(OUT / "fresh/estimates.npz"))
metrics = pd.read_csv(OUT / "fresh/metrics.tsv", sep="\t", low_memory=False)
samples = pd.read_csv(OUT / "fresh/samples.tsv", sep="\t")
selected = (
    metrics[metrics.method == decision["benchmark_label"]].groupby(["condition", "burden"]).head(1)
)
indices = dict(zip(samples["sample"], samples.index))
for index, row in selected.iterrows():
    i = indices[row["sample"]]
    info = {}
    weights, _, _ = decompose_mutational_profile_counts(
        inputs["counts"][i],
        (inputs["reference"], inputs["names"].tolist()),
        others_threshold=0,
        diagnostics=info,
    )
    np.testing.assert_allclose(
        weights, np.r_[estimates["known"][index], estimates["null"][index]], atol=1e-10, rtol=0
    )
    assert info["background_selection"] == DEFAULT_BACKGROUND_SELECTION
    assert info["fallback"] is None
result = {
    "method": DEFAULT_METHOD,
    "background_selection": DEFAULT_BACKGROUND_SELECTION,
    "matched_benchmark_fits": len(selected),
    "all_condition_burden_cells": True,
}
(OUT / "default-verification.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
