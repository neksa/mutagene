"""Test a frozen gradual search transition on an independent scenario panel."""

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1]))
sys.path.insert(0, str(ROOT.parent / "low-burden-2026"))
sys.path.insert(0, str(ROOT.parent / "signature-benchmark-2026"))

import frozen_sparse as sparse  # noqa: E402
from panel import generate_panel  # noqa: E402
from run_benchmark import score, sha256  # noqa: E402

from mutagene.signatures import get_dummy_signatures_lists  # noqa: E402


def blend(counts, current, forward):
    burden = counts.sum()
    if burden <= 1000 or burden >= 2000:
        return current.copy()
    fraction = (burden - 1000) / 1000
    return (1 - fraction) * forward + fraction * current


def evaluate(task):
    counts, truth, matrix, names, generating, meta = task
    current = sparse.fit_sparse(matrix, counts, len(names))
    forward = (
        current if counts.sum() <= 1000 else sparse._fit_forward(matrix, counts, len(names), None)
    )
    averaged = blend(counts, current, forward)
    if counts.sum() <= 1000 or counts.sum() >= 2000:
        np.testing.assert_array_equal(averaged, current)
    rows = []
    for method, w in [
        ("current", current),
        ("transition-average", averaged),
        ("forward-all", forward),
    ]:
        metrics = score(
            w[: len(names)],
            w[len(names) :],
            truth,
            counts,
            names,
            generating,
            matrix[:, : len(names)],
            meta["null_truth"],
        )
        rows.append(dict(meta, method=method, **metrics))
    return rows, [current, averaged, forward]


def main():
    out = ROOT / "transition"
    out.mkdir()
    spec = dict(
        seed=202610299,
        samples_per_cancer=2,
        burdens=[1000, 1001, 1250, 1500, 1750, 1999, 2000, 50000],
        formula="For 1000 < N < 2000, use (1-alpha)*forward + alpha*current, alpha=(N-1000)/1000. Otherwise return current.",
        gate="Both mean errors must improve within the transition range; preserve exact estimates at endpoint controls. Report condition regressions and paired intervals.",
        caveat="The blended estimate is not a constrained MLE refit or a BIC-selected subset.",
        sources={str(p): sha256(p) for p in [Path(__file__), ROOT / "frozen_sparse.py"]},
    )
    (out / "specification.json").write_text(json.dumps(spec, indent=2) + "\n")
    data, meta = generate_panel(spec["seed"], 2, spec["burdens"])
    meta.to_csv(out / "samples.tsv", sep="\t", index=False)
    np.savez_compressed(out / "inputs.npz", **data)
    matrix = np.column_stack(
        [data["reference"], np.array([v for _, v in get_dummy_signatures_lists()]).T]
    )
    tasks = [
        (
            data["counts"][i],
            data["truth"][i],
            matrix,
            data["names"].tolist(),
            data["generating_names"].tolist(),
            meta.iloc[i].to_dict(),
        )
        for i in range(len(meta))
    ]
    rows = []
    weights = []
    with ProcessPoolExecutor(max_workers=4) as pool:
        for done, (group, w) in enumerate(pool.map(evaluate, tasks), 1):
            rows.extend(group)
            weights.extend(w)
            if done % 16 == 0:
                pd.DataFrame(rows).to_csv(out / "metrics.tsv", sep="\t", index=False)
                print(f"Completed {done}/{len(tasks)} transition profiles", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "metrics.tsv", sep="\t", index=False)
    np.savez_compressed(out / "estimates.npz", weights=np.array(weights))
    frame.groupby(["burden", "method"])[["aggregate_null_tae", "paper_tae"]].mean().to_csv(
        out / "summary.tsv", sep="\t"
    )
    print("Completed", len(frame), "transition estimates", flush=True)


if __name__ == "__main__":
    main()
