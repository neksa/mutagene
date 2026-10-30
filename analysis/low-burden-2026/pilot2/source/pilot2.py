"""Screen candidates on previously examined profiles; reserve fresh draws for validation."""

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(ROOT.parent / "signature-benchmark-2026"))

from candidates import ROUND2 as CONFIGS, fit_candidate  # noqa: E402
from run_benchmark import score, sha256  # noqa: E402

from mutagene.signatures import get_dummy_signatures_lists  # noqa: E402


def evaluate(task):
    index, matrix, counts, truth, names, generating, metadata = task
    rows = []
    for config in CONFIGS:
        start = time.perf_counter()
        weights = fit_candidate(matrix, counts, len(names), config)
        metrics = score(
            weights[: len(names)],
            weights[len(names) :],
            truth,
            counts,
            names,
            generating,
            matrix[:, : len(names)],
            metadata["null_truth"],
        )
        rows.append(
            dict(metadata, method=config[0], seconds=time.perf_counter() - start, **metrics)
        )
    return rows


def main():
    out = ROOT / "pilot2"
    out.mkdir()
    source = ROOT.parent / "signature-benchmark-2026/background-results/paired"
    data = dict(np.load(source / "inputs.npz"))
    metadata = pd.read_csv(source / "samples.tsv", sep="\t")
    indices = metadata.index[(metadata.burden == 100) & (metadata.sample_index < 2)]
    dummy = np.array([v for _, v in get_dummy_signatures_lists()]).T
    matrix = np.column_stack([data["reference"], dummy])
    (out / "specification.json").write_text(
        json.dumps(
            {
                "purpose": "Exploratory screening only; these profiles already inform earlier work.",
                "profiles": len(indices),
                "configs": CONFIGS,
                "input_sha256": sha256(source / "inputs.npz"),
                "candidate_sha256": sha256(ROOT / "candidates.py"),
                "runner_sha256": sha256(__file__),
                "selection": "At 100 mutations, sample_index < 2; all cancers and conditions.",
            },
            indent=2,
        )
        + "\n"
    )
    tasks = [
        (
            i,
            matrix,
            data["counts"][i],
            data["truth"][i],
            data["names"].tolist(),
            data["generating_names"].tolist(),
            metadata.iloc[i].to_dict(),
        )
        for i in indices
    ]
    rows = []
    with ProcessPoolExecutor(max_workers=4) as pool:
        for done, group in enumerate(pool.map(evaluate, tasks), 1):
            rows.extend(group)
            pd.DataFrame(rows).to_csv(out / "metrics.tsv", sep="\t", index=False)
            print(f"Completed {done}/{len(tasks)} profiles", flush=True)
    frame = pd.DataFrame(rows)
    summary = frame.groupby("method")[
        ["aggregate_null_tae", "paper_tae", "precision", "recall", "seconds"]
    ].mean()
    summary.to_csv(out / "summary.tsv", sep="\t")
    print(summary.sort_values("aggregate_null_tae").to_string(), flush=True)


if __name__ == "__main__":
    main()
