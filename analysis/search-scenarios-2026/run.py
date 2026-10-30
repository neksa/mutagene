"""Compare search refinements across benchmark scenarios without changing the default."""

import json
import sys
import time
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
from mutagene.signatures.likelihood import fit_likelihood  # noqa: E402


def swaps(matrix, counts, known_count, weights, info):
    selected = tuple(
        sorted(
            info["selected_known_indices"]
            + [known_count + i for i in info["selected_background_indices"]]
        )
    )
    current = (weights, info["criterion_value"])
    cache = {selected: current}
    for step in range(3):
        best = current
        chosen = selected
        for remove in selected:
            for add in range(matrix.shape[1]):
                if add in selected:
                    continue
                subset = tuple(sorted([i for i in selected if i != remove] + [add]))
                if subset not in cache:
                    if (matrix[counts > 0][:, subset].sum(axis=1) == 0).any():
                        continue
                    diagnostic = {}
                    w = fit_likelihood(matrix[:, subset], counts, diagnostic)
                    full = np.zeros(matrix.shape[1])
                    full[list(subset)] = w
                    value = sparse.information_criterion(
                        diagnostic["nll_per_mutation"] * counts.sum(),
                        len(subset) - 1,
                        counts.sum(),
                        "bic",
                    )
                    cache[subset] = (full, value)
                if cache[subset][1] < best[1] - 1e-7:
                    best, chosen = cache[subset], subset
        if chosen == selected:
            break
        selected, current = chosen, best
    return current[0], len(cache)


def evaluate(task):
    counts, truth, matrix, names, generating, metadata = task
    start = time.perf_counter()
    info = {}
    baseline = sparse.fit_sparse(matrix, counts, len(names), diagnostics=info)
    full_info = {}
    forward_info = {}
    if counts.sum() <= 1000:
        forward, forward_info = baseline, info
        full = sparse.fit_sparse(
            matrix, counts, len(names), diagnostics=full_info, search_start="full"
        )
    else:
        full, full_info = baseline, info
        forward = sparse._fit_forward(matrix, counts, len(names), forward_info)
    values = np.array([full_info["criterion_value"], forward_info["criterion_value"]])
    best = [full, forward][int(np.argmin(values))]
    distinct = not np.allclose(full, forward, atol=1e-10, rtol=0)
    probabilities = np.exp(-0.5 * (values - values.min()))
    probabilities /= probabilities.sum()
    average = probabilities[0] * full + probabilities[1] * forward
    swapped, refits = swaps(matrix, counts, len(names), baseline, info)
    fixed = sparse.fit_sparse(matrix, counts, len(names), background_selection="fixed")
    outputs = [
        ("current", baseline),
        ("fixed-BIC", fixed),
        ("forward-all", forward),
        ("best-start-BIC", best),
        ("BIC-weighted-average", average),
        ("swap-BIC", swapped),
    ]
    rows = []
    for method, weights in outputs:
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
            dict(
                metadata,
                method=method,
                distinct_starts=distinct,
                swap_refits=refits,
                full_bic=values[0],
                forward_bic=values[1],
                **metrics,
            )
        )
    return rows, [w for _, w in outputs], time.perf_counter() - start


def main():
    out = ROOT / "screen"
    out.mkdir()
    spec = dict(
        seed=202610296,
        samples_per_cancer=2,
        burdens=[100, 250, 500, 999, 1000, 1001, 2000, 50000],
        sensitivity_seed=202610297,
        sensitivity_samples_per_cancer=1,
        sensitivity_burdens=[50, 100, 1000],
        methods=[
            "current",
            "fixed-BIC",
            "forward-all",
            "best-start-BIC",
            "BIC-weighted-average",
            "swap-BIC",
        ],
        swap_policy="At most three greedy one-for-one replacement rounds; retain model size.",
        averaging="Weight two distinct search results by exp(-delta_BIC/2); this is not a posterior over all models.",
        decision="Screening only. Require improvement in both low-burden metrics without mean high-burden regression before independent confirmation. Retain all scenarios and failures.",
        sources={
            str(p): sha256(p)
            for p in [
                Path(__file__),
                ROOT / "frozen_sparse.py",
                ROOT.parent / "low-burden-2026/panel.py",
                ROOT.parents[1] / "mutagene/signatures/likelihood.py",
            ]
        },
    )
    (out / "specification.json").write_text(json.dumps(spec, indent=2) + "\n")
    primary, meta = generate_panel(spec["seed"], 2, spec["burdens"])
    stress, sm = generate_panel(spec["sensitivity_seed"], 1, spec["sensitivity_burdens"], False)
    meta["panel"] = "primary"
    sm["panel"] = "no-truth-cutoff"
    metadata = pd.concat([meta, sm], ignore_index=True)
    data = dict(primary)
    data["counts"] = np.concatenate([primary["counts"], stress["counts"]])
    data["truth"] = np.concatenate([primary["truth"], stress["truth"]])
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
    rows = []
    estimates = []
    times = []
    with ProcessPoolExecutor(max_workers=4) as pool:
        for done, (group, weights, seconds) in enumerate(pool.map(evaluate, tasks), 1):
            rows.extend(group)
            estimates.extend(weights)
            times.append(seconds)
            if done % 16 == 0:
                pd.DataFrame(rows).to_csv(out / "metrics.tsv", sep="\t", index=False)
                print(f"Completed {done}/{len(tasks)} profiles", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "metrics.tsv", sep="\t", index=False)
    np.savez_compressed(out / "estimates.npz", weights=np.array(estimates))
    frame.groupby(["panel", "burden", "method"])[["aggregate_null_tae", "paper_tae"]].mean().to_csv(
        out / "summary.tsv", sep="\t"
    )
    (out / "completion.json").write_text(
        json.dumps(
            dict(
                profiles=len(metadata),
                fits=len(frame),
                errors=0,
                mean_profile_seconds=float(np.mean(times)),
            ),
            indent=2,
        )
        + "\n"
    )
    print("Complete", len(frame), "fits", flush=True)


if __name__ == "__main__":
    main()
