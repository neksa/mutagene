"""Isolate the search-policy change after adding one mutation to a 1,000-mutation profile."""

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1]))
sys.path.insert(0, str(ROOT.parent / "signature-benchmark-2026"))

import frozen_sparse as sparse  # noqa: E402
from run_benchmark import score  # noqa: E402

from mutagene.signatures import get_dummy_signatures_lists  # noqa: E402


def evaluate(task):
    index, counts, truth, matrix, names, generating, meta = task
    rng = np.random.default_rng(np.random.SeedSequence([202610298, index]))
    channel = int(rng.choice(len(counts), p=counts / counts.sum()))
    incremented = counts.copy()
    incremented[channel] += 1
    before = sparse.fit_sparse(matrix, counts, len(names))
    after = sparse.fit_sparse(matrix, incremented, len(names))
    forward = sparse._fit_forward(matrix, incremented, len(names), None)

    def metric(w, c):
        return score(
            w[: len(names)],
            w[len(names) :],
            truth,
            c,
            names,
            generating,
            matrix[:, : len(names)],
            meta["null_truth"],
        )["aggregate_null_tae"]

    return dict(
        sample=meta["sample"],
        cancer=meta["cancer"],
        condition=meta["condition"],
        added_channel=channel,
        current_l1_change=float(np.abs(after - before).sum()),
        forward_l1_change=float(np.abs(forward - before).sum()),
        same_counts_search_l1=float(np.abs(after - forward).sum()),
        error_1000=metric(before, counts),
        error_1001=metric(after, incremented),
        forward_error_1001=metric(forward, incremented),
    )


def main():
    out = ROOT / "boundary"
    out.mkdir()
    data = dict(np.load(ROOT / "screen/inputs.npz"))
    metadata = pd.read_csv(ROOT / "screen/samples.tsv", sep="\t")
    indices = metadata.index[(metadata.panel == "primary") & (metadata.burden == 1000)]
    matrix = np.column_stack(
        [data["reference"], np.array([v for _, v in get_dummy_signatures_lists()]).T]
    )
    tasks = [
        (
            i,
            data["counts"][i],
            data["truth"][i],
            matrix,
            data["names"].tolist(),
            data["generating_names"].tolist(),
            metadata.iloc[i].to_dict(),
        )
        for i in indices
    ]
    with ProcessPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(evaluate, tasks))
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "metrics.tsv", sep="\t", index=False)
    summary = {
        name: float(frame[name].mean())
        for name in [
            "current_l1_change",
            "forward_l1_change",
            "same_counts_search_l1",
            "error_1000",
            "error_1001",
            "forward_error_1001",
        ]
    }
    summary.update(
        profiles=len(frame),
        max_current_l1_change=float(frame.current_l1_change.max()),
        seed=202610298,
        interpretation="Add one mutation sampled from each observed profile; compare both searches on the identical 1001-mutation counts.",
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
