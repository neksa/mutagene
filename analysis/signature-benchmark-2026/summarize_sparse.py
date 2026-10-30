"""Summarize the frozen sparse evaluation with visible reconstruction and attribution metrics."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import vl_convert as vlc

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "sparse-results"
FIGURES = OUT / "figures"
FIGURES.mkdir(exist_ok=True)
metrics = pd.read_csv(OUT / "evaluate/metrics.tsv", sep="\t", low_memory=False)
if len(metrics) != 768 * 27:
    raise ValueError("Evaluation must complete before report generation")
selection = json.loads((OUT / "selection.json").read_text())
selected = selection["selected"]
metrics["method_label"] = metrics.method + np.where(
    metrics.background, " + background", " / no background"
)
metrics.loc[metrics.method == "SigProfilerAssignment", "method_label"] = "SigProfilerAssignment"
metrics["family"] = np.select(
    [
        metrics.method == selected,
        metrics.method.str.startswith("sparse-"),
        metrics.method == "SigProfilerAssignment",
    ],
    ["Selected sparse criterion", "Other sparse criterion", "SigProfilerAssignment"],
    default="Existing objective",
)
metrics["failed"] = metrics.error.notna()
metrics["fallback_used"] = metrics.fallback.notna()
summary = (
    metrics.groupby(["condition", "burden", "method", "background", "method_label", "family"])
    .mean(numeric_only=True)
    .reset_index()
)
summary.to_csv(OUT / "evaluation-summary.tsv", sep="\t", index=False)
records = json.loads(summary.to_json(orient="records"))
colors = {
    "field": "family",
    "type": "nominal",
    "title": None,
    "scale": {
        "domain": [
            "Selected sparse criterion",
            "Other sparse criterion",
            "SigProfilerAssignment",
            "Existing objective",
        ],
        "range": ["#00796b", "#e69f00", "#b32879", "#7c8792"],
    },
}
common = {
    "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
    "data": {"values": records},
    "params": [
        {
            "name": "condition_choice",
            "value": "clean",
            "bind": {
                "input": "select",
                "options": ["clean", "out_of_reference_20", "uniform_20", "substitution_noise_20"],
                "name": "Condition: ",
            },
        },
        {
            "name": "burden_choice",
            "value": 2000,
            "bind": {"input": "select", "options": [100, 2000, 50000], "name": "Mutations: "},
        },
    ],
    "transform": [
        {"filter": "datum.condition === condition_choice && datum.burden === burden_choice"}
    ],
    "config": {
        "font": "Arial",
        "view": {"stroke": None},
        "legend": {"orient": "bottom", "columns": 2},
        "axis": {"labelFontSize": 11, "titleFontSize": 12},
    },
}
row = {
    "field": "method_label",
    "type": "nominal",
    "sort": {"field": "aggregate_null_tae", "order": "ascending"},
    "title": None,
    "axis": {"labelLimit": 300},
}
tooltip = [
    {"field": f, "format": ".4f"}
    for f in [
        "aggregate_null_tae",
        "cosine",
        "active_known",
        "precision",
        "recall",
        "null_weight",
        "fallback_used",
        "failed",
    ]
]
chart = {
    **common,
    "hconcat": [
        {
            "width": 270,
            "height": 560,
            "mark": {"type": "bar"},
            "encoding": {
                "y": row,
                "x": {
                    "field": "aggregate_null_tae",
                    "type": "quantitative",
                    "title": "Attribution error (lower is better)",
                    "scale": {"domain": [0, 1]},
                },
                "color": colors,
                "tooltip": tooltip,
            },
        },
        {
            "width": 240,
            "height": 560,
            "mark": {"type": "point", "filled": True, "size": 65},
            "encoding": {
                "y": {**row, "axis": None},
                "x": {
                    "field": "cosine",
                    "type": "quantitative",
                    "title": "Reconstruction cosine (higher is better)",
                    "scale": {"zero": False},
                },
                "color": colors,
                "tooltip": tooltip,
            },
        },
    ],
    "resolve": {"scale": {"y": "shared"}},
}
(FIGURES / "all-methods.vl.json").write_text(json.dumps(chart, indent=2) + "\n")
(FIGURES / "all-methods.svg").write_text(vlc.vegalite_to_svg(chart))
(FIGURES / "all-methods.png").write_bytes(vlc.vegalite_to_png(chart, scale=1.5))
(FIGURES / "interactive.html").write_text(vlc.vegalite_to_html(chart, bundle=True))

rng = np.random.default_rng(202610233)
paired = []
for (condition, burden), frame in metrics.groupby(["condition", "burden"]):
    for background in [True, False]:
        candidate = frame[(frame.method == selected) & (frame.background == background)].set_index(
            ["cancer", "sample"]
        )
        for comparator in ["MLE", "SigProfilerAssignment"]:
            baseline = frame[
                (frame.method == comparator)
                & (frame.background == (background if comparator == "MLE" else False))
            ].set_index(["cancer", "sample"])
            delta = candidate.aggregate_null_tae - baseline.aggregate_null_tae
            if delta.isna().any():
                raise ValueError("Failed paired fits require explicit handling")
            bootstrap = np.zeros(2000)
            for _, group in delta.groupby(level="cancer"):
                values = group.to_numpy()
                bootstrap += (
                    values[rng.integers(len(values), size=(2000, len(values)))].mean(axis=1) / 8
                )
            paired.append(
                {
                    "condition": condition,
                    "burden": burden,
                    "background": background,
                    "comparator": comparator,
                    "sparse_minus_comparator": delta.mean(),
                    "ci_low": np.quantile(bootstrap, 0.025),
                    "ci_high": np.quantile(bootstrap, 0.975),
                }
            )
pd.DataFrame(paired).to_csv(OUT / "paired-differences.tsv", sep="\t", index=False)

clean = summary[
    (summary.condition == "clean")
    & (summary.method.isin(["MLE", selected, "SigProfilerAssignment"]))
]
lines = [
    "| Mutations | Method | Attribution error | Reconstruction cosine | Active signatures | Precision | Recall |",
    "|---:|---|---:|---:|---:|---:|---:|",
]
for _, r in clean.sort_values(["burden", "method_label"]).iterrows():
    lines.append(
        f"| {r.burden} | {r.method_label} | {r.aggregate_null_tae:.3f} | {r.cosine:.6f} | {r.active_known:.2f} | {r.precision:.3f} | {r.recall:.3f} |"
    )
(OUT / "clean-comparison.md").write_text("\n".join(lines) + "\n")
print(
    clean[
        [
            "burden",
            "method_label",
            "aggregate_null_tae",
            "active_known",
            "precision",
            "recall",
            "cosine",
        ]
    ].to_string(index=False)
)
print("Failed fits:", int(metrics.failed.sum()))
