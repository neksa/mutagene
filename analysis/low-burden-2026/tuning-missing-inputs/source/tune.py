"""Compare subset-search initialization across mutation burdens on fresh tuning draws."""

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1]))
sys.path.insert(0, str(ROOT.parent / "signature-benchmark-2026"))

from candidates import fit_candidate  # noqa: E402
from panel import generate_panel  # noqa: E402
from run_benchmark import score, sha256  # noqa: E402

from mutagene.signatures import get_dummy_signatures_lists  # noqa: E402
from mutagene.signatures.likelihood import negative_log_likelihood  # noqa: E402

CONFIGS = [
    ("joint-BIC", "bic", 0.0, True, False),
    ("forward-BIC", "forward", 0.0, True, False),
    ("forward-EBIC-0.5", "forward", 0.5, True, False),
]


def evaluate(task):
    counts, truth, matrix, names, generating, metadata = task
    rows, estimates = [], []
    for config in CONFIGS:
        begin = time.perf_counter()
        weights = fit_candidate(matrix, counts, len(names), config)
        elapsed = time.perf_counter() - begin
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
        criterion = 2 * negative_log_likelihood(weights, matrix, counts) + (
            np.count_nonzero(weights > 1e-7) - 1
        ) * np.log(counts.sum())
        rows.append(
            dict(metadata, method=config[0], seconds=elapsed, bic=float(criterion), **metrics)
        )
        estimates.append(weights)
    return rows, estimates


def main():
    out = ROOT / "tuning"
    out.mkdir()
    seed, samples, burdens = 202610291, 4, [100, 250, 500, 1000, 2000, 50000]
    specification = dict(
        seed=seed,
        samples_per_cancer=samples,
        burdens=burdens,
        configs=CONFIGS,
        purpose="Select search initialization and any burden boundary before fresh validation.",
        source_sha256={
            name: sha256(ROOT / name) for name in ["tune.py", "panel.py", "candidates.py"]
        },
        selection="Require improvement in both low-burden errors; report every condition and high-burden regression.",
    )
    (out / "specification.json").write_text(json.dumps(specification, indent=2) + "\n")
    data, metadata = generate_panel(seed, samples, burdens)
    metadata.to_csv(out / "samples.tsv", sep="\t", index=False)
    np.savez_compressed(out / "inputs.npz", **data)
    dummy = np.array([v for _, v in get_dummy_signatures_lists()]).T
    matrix = np.column_stack([data["reference"], dummy])
    tasks = [
        (
            data["counts"][i],
            data["truth"][i],
            matrix,
            data["names"].tolist(),
            data["generating_names"].tolist(),
            metadata.iloc[i].to_dict(),
        )
        for i in range(len(metadata))
    ]
    rows, estimates = [], []
    with ProcessPoolExecutor(max_workers=4) as pool:
        for done, (group, weights) in enumerate(pool.map(evaluate, tasks), 1):
            rows.extend(group)
            estimates.extend(weights)
            if done % 16 == 0:
                pd.DataFrame(rows).to_csv(out / "metrics.tsv", sep="\t", index=False)
                print(f"Completed {done}/{len(tasks)} profiles", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "metrics.tsv", sep="\t", index=False)
    np.savez_compressed(out / "estimates.npz", weights=np.array(estimates))
    summary = frame.groupby(["burden", "method"])[
        ["aggregate_null_tae", "paper_tae", "bic", "seconds"]
    ].mean()
    summary.to_csv(out / "summary.tsv", sep="\t")
    print(summary.to_string(), flush=True)


if __name__ == "__main__":
    main()
