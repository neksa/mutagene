"""Assess the independent transition experiment without changing its frozen rule."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def main():
    m = pd.read_csv(ROOT / "transition/metrics.tsv", sep="\t")
    assert len(m) == 1536
    baseline = m[m.method == "current"].set_index("sample")
    new = m[m.method == "transition-average"].copy()
    for metric in ["aggregate_null_tae", "paper_tae"]:
        new[metric] -= new["sample"].map(baseline[metric])
    new.groupby(["burden", "condition"])[["aggregate_null_tae", "paper_tae"]].mean().to_csv(
        ROOT / "results/transition-cell-differences.tsv", sep="\t"
    )
    rng = np.random.default_rng(202610300)
    rows = []
    for burden, g in new.groupby("burden"):
        for metric in ["aggregate_null_tae", "paper_tae"]:
            differences = g.groupby("cancer")[metric].mean().to_numpy()
            draws = rng.choice(differences, (4000, len(differences)), replace=True).mean(axis=1)
            lo, hi = np.quantile(draws, [0.025, 0.975])
            rows.append(
                dict(
                    burden=int(burden),
                    metric=metric,
                    difference=float(g[metric].mean()),
                    low=float(lo),
                    high=float(hi),
                )
            )
    pd.DataFrame(rows).to_csv(ROOT / "results/transition-intervals.tsv", sep="\t", index=False)
    meta = pd.read_csv(ROOT / "transition/samples.tsv", sep="\t")
    weights = np.load(ROOT / "transition/estimates.npz")["weights"].reshape(len(meta), 3, -1)
    endpoints = (meta.burden <= 1000) | (meta.burden >= 2000)
    exact = np.array_equal(weights[endpoints, 0], weights[endpoints, 1])
    middle = new[(new.burden > 1000) & (new.burden < 2000)]
    result = dict(
        profiles=len(meta),
        exact_endpoint_profiles=int(endpoints.sum()),
        exact_endpoint_match=exact,
        transition_total_error_change=float(middle.aggregate_null_tae.mean()),
        transition_thresholded_error_change=float(middle.paper_tae.mean()),
    )
    (ROOT / "results/transition-gates.json").write_text(json.dumps(result, indent=2) + "\n")
    summary = (
        m[m.method != "forward-all"]
        .groupby(["burden", "method"])[["aggregate_null_tae", "paper_tae"]]
        .mean()
    )
    table = [
        "| Mutations | Current total error | Transition total error | Current thresholded error | Transition thresholded error |",
        "|---|---:|---:|---:|---:|",
    ]
    for n in sorted(m.burden.unique()):
        a = summary.loc[(n, "current")]
        b = summary.loc[(n, "transition-average")]
        table.append(
            f"| {n:,} | {a.iloc[0]:.5f} | {b.iloc[0]:.5f} | {a.iloc[1]:.5f} | {b.iloc[1]:.5f} |"
        )
    (ROOT / "results/transition-tables.md").write_text("\n".join(table) + "\n")
    print(json.dumps(result, indent=2))
    print(pd.DataFrame(rows).to_string(index=False))
    print("\n".join(table))


if __name__ == "__main__":
    main()
