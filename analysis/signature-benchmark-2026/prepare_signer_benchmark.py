"""Freeze a balanced, matched signeR extension before evaluating its attribution accuracy."""

import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from run_benchmark import sha256

from mutagene.signatures import get_dummy_signatures_lists

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "signer-results"
SOURCE = ROOT / "background-results/fresh"
if OUT.exists():
    raise ValueError("Preserve existing evidence; choose an isolated checkout for a repeat run")
OUT.mkdir()
(OUT / "fits").mkdir()
inputs = dict(np.load(SOURCE / "inputs.npz"))
metadata = pd.read_csv(SOURCE / "samples.tsv", sep="\t")
selected = metadata[metadata.sample_index == 0].copy()
assert len(selected) == 96
assert selected.groupby(["cancer", "condition", "burden"]).size().eq(1).all()
selected["input_index"] = selected.index
selected.to_csv(OUT / "samples.tsv", sep="\t", index=False)
reference = pd.read_csv(SOURCE / "reference.tsv", sep="\t", index_col=0)
reference.to_csv(OUT / "reference.tsv", sep="\t")
counts = inputs["counts"][selected.index]
pd.DataFrame(counts.T, index=reference.index, columns=selected["sample"]).to_csv(
    OUT / "catalogs.tsv", sep="\t", index_label="MutationType"
)
background = get_dummy_signatures_lists()
pd.DataFrame(
    np.array([v for _, v in background]).T,
    index=reference.index,
    columns=["background_" + n for n, _ in background],
).to_csv(OUT / "background.tsv", sep="\t")
tasks = []
for row in selected.itertuples():
    for enabled in (False, True):
        for chain in ([1, 2] if row.cancer == "Liver-HCC" else [1]):
            task_id = len(tasks) + 1
            tasks.append(
                dict(
                    task_id=task_id,
                    sample=row.sample,
                    background=enabled,
                    chain=chain,
                    seed=202610280 + task_id,
                    diagnostic=row.cancer == "Liver-HCC",
                )
            )
pd.DataFrame(tasks).to_csv(OUT / "tasks.tsv", sep="\t", index=False)
specification = {
    "purpose": "Compare native fixed-catalog signeR inference on the same profiles as all existing comparators.",
    "profiles": 96,
    "primary_fits": 192,
    "independent_repeat_fits": 24,
    "selection": "sample_index == 0 in every cancer/condition/burden cell; no selection by fit accuracy",
    "phase": "subset of background-results/fresh, seed 202610273",
    "methods": ["signeR", "signeR+background"],
    "package": "signeR",
    "version": "2.14.0",
    "upstream_commit": subprocess.check_output(
        ["git", "-C", "/tmp/mutagene-signer-source", "rev-parse", "HEAD"], text=True
    ).strip(),
    "settings": {
        "fixedP": True,
        "main_burn": 10000,
        "main_eval": 2000,
        "estimate_hyper": False,
        "opportunities": "ones; benchmark reference already represents channel probabilities",
        "parallelization": "sequential within fit; eight independent workers",
    },
    "prior": "Native 2.14.0 defaults; no hyperparameter tuning on benchmark truth.",
    "estimator": "Native Median_exp, then divide by its total to compare relative attribution. Preserve raw median counts and raw total/burden.",
    "background_variant": "Append the same six context-independent components. This is a configured reference extension, not a native signeR null-test feature.",
    "diagnostics": "Single-chain split Rhat and ESS for every fit; independent second chains for both variants in all 12 Liver-HCC condition/burden cells. Do not equate numerical completion with convergence.",
    "ranking": "Primary: aggregate_null_tae. Also report paper_tae, background mass, cosine, active known count, runtime, raw-scale sensitivity, and chain stability.",
    "input_sha256": {
        name: sha256(SOURCE / name)
        for name in ["inputs.npz", "samples.tsv", "reference.tsv", "metrics.tsv", "estimates.npz"]
    },
    "adapter_sha256": sha256(ROOT / "run_signer.R"),
    "preparation_sha256": sha256(__file__),
    "limitations": [
        "A 96-profile stratified subset, not the full 768-profile panel.",
        "Fixed-catalog attribution does not benchmark de novo discovery.",
        "Single-profile fits use native fixed hyperparameters; they do not borrow information across profiles.",
        "Existing synthetic opportunities do not validate opportunity correction on real genomic targets.",
    ],
}
(OUT / "specification.json").write_text(json.dumps(specification, indent=2) + "\n")
print(f"Prepared {len(selected)} profiles and {len(tasks)} fits")
