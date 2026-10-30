"""Report scenario-level gains and regressions against the current default."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import vl_convert as vlc

ROOT = Path(__file__).resolve().parent
METRICS = ["aggregate_null_tae", "paper_tae"]


def main():
    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    m = pd.read_csv(ROOT / "screen/metrics.tsv", sep="\t")
    assert len(m) == 3648 and not m.duplicated(["sample", "method"]).any()
    summary = m.groupby(["panel", "burden", "method"])[METRICS].mean().reset_index()
    summary.to_csv(out / "summary.tsv", sep="\t", index=False)
    m.groupby(["panel", "burden", "condition", "method"])[METRICS].mean().reset_index().to_csv(
        out / "conditions.tsv", sep="\t", index=False
    )
    m.groupby(["panel", "burden", "cancer", "method"])[METRICS].mean().reset_index().to_csv(
        out / "cancers.tsv", sep="\t", index=False
    )
    baseline = m[m.method == "current"].set_index("sample")
    diffs = m[m.method != "current"].copy()
    for metric in METRICS:
        diffs[metric] = diffs[metric] - diffs["sample"].map(baseline[metric])
    diffs.to_csv(out / "paired-differences.tsv", sep="\t", index=False)
    cells = diffs.groupby(["panel", "burden", "condition", "method"])[METRICS].mean().reset_index()
    cells.to_csv(out / "cell-differences.tsv", sep="\t", index=False)
    rng = np.random.default_rng(202610298)
    comparisons = []
    for (panel, method), group in diffs.groupby(["panel", "method"]):
        for label, subset in [
            ("low", group[group.burden <= 1000]),
            ("high", group[group.burden > 1000]),
        ]:
            if subset.empty:
                continue
            for metric in METRICS:
                by_cancer = subset.groupby("cancer")[metric].mean().to_numpy()
                draws = rng.choice(by_cancer, (4000, len(by_cancer)), replace=True).mean(axis=1)
                lo, hi = np.quantile(draws, [0.025, 0.975])
                comparisons.append(
                    dict(
                        panel=panel,
                        method=method,
                        range=label,
                        metric=metric,
                        difference=subset[metric].mean(),
                        low=lo,
                        high=hi,
                    )
                )
    paired = pd.DataFrame(comparisons)
    paired.to_csv(out / "intervals.tsv", sep="\t", index=False)
    verdicts = {}
    for method, g in paired[paired.panel == "primary"].groupby("method"):
        low = g[g["range"] == "low"]
        high = g[g["range"] == "high"]
        verdicts[method] = dict(
            low_mean_improves_both=bool((low.difference < 0).all()),
            low_intervals_below_zero=bool((low.high < 0).all()),
            high_mean_no_regression=bool((high.difference <= 1e-12).all()),
        )
    (out / "screening-gates.json").write_text(json.dumps(verdicts, indent=2) + "\n")
    text = [
        "| Method | Low-burden total error | Low-burden thresholded error | High-burden total error | High-burden thresholded error |",
        "|---|---:|---:|---:|---:|",
    ]
    primary = m[m.panel == "primary"]
    for method, g in primary.groupby("method"):
        low = g[g.burden <= 1000][METRICS].mean()
        high = g[g.burden > 1000][METRICS].mean()
        text.append(
            f"| {method} | {low.iloc[0]:.5f} | {low.iloc[1]:.5f} | {high.iloc[0]:.5f} | {high.iloc[1]:.5f} |"
        )
    (out / "tables.md").write_text("\n".join(text) + "\n")
    chart = cells[(cells.panel == "primary")].copy()
    spec = {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "data": {"values": chart.to_dict("records")},
        "facet": {"row": {"field": "condition", "type": "nominal", "title": None}},
        "spec": {
            "width": 680,
            "height": 130,
            "mark": "rect",
            "encoding": {
                "x": {
                    "field": "burden",
                    "type": "ordinal",
                    "title": "Mutation burden",
                    "sort": [100, 250, 500, 999, 1000, 1001, 2000, 50000],
                },
                "y": {"field": "method", "type": "nominal", "title": None},
                "color": {
                    "field": "aggregate_null_tae",
                    "type": "quantitative",
                    "scale": {"scheme": "redblue", "domainMid": 0, "reverse": True},
                    "title": "Error change",
                },
                "tooltip": [
                    {"field": x}
                    for x in ["condition", "burden", "method", "aggregate_null_tae", "paper_tae"]
                ],
            },
        },
        "config": {"font": "Arial", "view": {"stroke": None}},
    }
    (out / "comparison.vl.json").write_text(json.dumps(spec, indent=2) + "\n")
    (out / "comparison.svg").write_text(vlc.vegalite_to_svg(spec))
    (out / "comparison.png").write_bytes(vlc.vegalite_to_png(spec, scale=2))
    (out / "interactive.html").write_text(vlc.vegalite_to_html(spec, bundle=True))
    print(json.dumps(verdicts, indent=2))
    print("\n".join(text))


if __name__ == "__main__":
    main()
