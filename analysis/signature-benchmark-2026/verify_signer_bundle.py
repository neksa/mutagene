"""Verify the finalized signeR bundle and retained historical artifact manifests."""

import hashlib
import json
from pathlib import Path

import pandas as pd
from signer_adapter import read_fit

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
OUT = ROOT / "signer-results"


def verify_manifest(path):
    records = json.loads(path.read_text())
    for name, record in records.items():
        artifact = REPO / name
        assert artifact.stat().st_size == record["bytes"], name
        assert hashlib.sha256(artifact.read_bytes()).hexdigest() == record["sha256"], name
    return len(records)


def main():
    tasks = pd.read_csv(OUT / "tasks.tsv", sep="\t")
    reference = pd.read_csv(OUT / "reference.tsv", sep="\t", index_col=0)
    background = pd.read_csv(OUT / "background.tsv", sep="\t", index_col=0)
    expected_ids = set(tasks.task_id)
    assert len(tasks) == len(expected_ids) == 216
    assert {int(path.stem) for path in (OUT / "fits").glob("*.json")} == expected_ids
    for task in tasks.to_dict("records"):
        prefix = OUT / "fits" / f"{task['task_id']:04d}"
        _, _, info = read_fit(prefix, task, reference.columns, background.columns)
        assert info["chain"] == task["chain"]
        assert info["seed"] == task["seed"]
        assert info["posterior_draws"] == 2000
        if task["diagnostic"]:
            assert prefix.with_suffix(".draws.rds").stat().st_size > 0
        assert prefix.with_suffix(".diagnostics.tsv").stat().st_size > 0

    pairs = pd.read_csv(OUT / "chain-stability.tsv", sep="\t")
    assert len(pairs) == 24
    assert set(pairs.repeat_task) == set(tasks.loc[tasks.chain == 2, "task_id"])
    for pair in pairs.itertuples():
        first = tasks.loc[tasks.task_id == pair.first_task].iloc[0]
        repeat = tasks.loc[tasks.task_id == pair.repeat_task].iloc[0]
        assert first.chain == 1 and repeat.chain == 2
        assert first["sample"] == repeat["sample"] == pair.sample
        assert first.background == repeat.background == pair.background
        assert (OUT / f"chain-diagnostics-{pair.first_task:04d}.tsv").stat().st_size > 0

    manifests = [
        ROOT / "artifact-manifest.json",
        *[
            ROOT / name / "artifact-manifest.json"
            for name in (
                "likelihood-audit",
                "sparse-results",
                "repaired-results",
                "background-results",
                "signer-results",
            )
        ],
    ]
    counts = {str(path.relative_to(REPO)): verify_manifest(path) for path in manifests}
    print(json.dumps({"fits": 216, "chain_pairs": 24, "manifest_entries": counts}, indent=2))


if __name__ == "__main__":
    main()
