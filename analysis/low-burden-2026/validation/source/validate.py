"""Evaluate the frozen burden policy on fresh profiles and a truth-cutoff sensitivity panel."""

import argparse
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


def evaluate(task):
    counts, truth, matrix, names, generating, metadata, selection = task
    baseline = ("joint-BIC", "bic", 0.0, True, False)
    candidate = selection["config"] if counts.sum() <= selection["maximum_burden"] else baseline
    rows, estimates = [], []
    for label, config in [("joint-BIC", baseline), ("adaptive-BIC", candidate)]:
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
        rows.append(
            dict(metadata, method=label, selected_search=config[1], seconds=elapsed, **metrics)
        )
        estimates.append(weights)
    if counts.sum() > selection["maximum_burden"]:
        np.testing.assert_array_equal(estimates[0], estimates[1])
    return rows, estimates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["validation", "sensitivity"])
    args = parser.parse_args()
    out = ROOT / args.phase
    out.mkdir()
    selection = json.loads((ROOT / "selection.json").read_text())
    if args.phase == "validation":
        seed, samples, burdens, cutoff = 202610292, 8, [100, 250, 500, 1000, 2000, 50000], True
    else:
        seed, samples, burdens, cutoff = 202610293, 4, [50, 100], False
    spec = dict(
        seed=seed,
        samples_per_cancer=samples,
        burdens=burdens,
        truth_cutoff=cutoff,
        selection_sha256=sha256(ROOT / "selection.json"),
        selection=selection,
        source_sha256={
            name: sha256(ROOT / name) for name in ["validate.py", "panel.py", "candidates.py"]
        },
        package_source_sha256={
            name: sha256(ROOT.parents[1] / "mutagene/signatures" / name)
            for name in ["sparse.py", "likelihood.py"]
        },
    )
    (out / "specification.json").write_text(json.dumps(spec, indent=2) + "\n")
    data, metadata = generate_panel(seed, samples, burdens, cutoff)
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
            selection,
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
        ["aggregate_null_tae", "paper_tae", "seconds"]
    ].mean()
    summary.to_csv(out / "summary.tsv", sep="\t")
    print(summary.to_string(), flush=True)


if __name__ == "__main__":
    main()
