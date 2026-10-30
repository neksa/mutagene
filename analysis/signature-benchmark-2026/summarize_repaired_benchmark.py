"""Report repaired attribution, reconstruction, and paired uncertainty without retuning."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import vl_convert as vlc

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "repaired-results"
FIGURES = OUT / "figures"
FIGURES.mkdir(exist_ok=True)
frames = []
for phase, expected in [("paired", 20736), ("fresh", 22272)]:
    frame = pd.read_csv(OUT / phase / "metrics.tsv", sep="\t", low_memory=False)
    if len(frame) != expected or frame.error.notna().any():
        raise ValueError(f"Incomplete or failed {phase} fits require investigation")
    if phase == "paired":
        archived = pd.read_csv(
            ROOT / "sparse-results/evaluate/metrics.tsv", sep="\t", low_memory=False
        )
        archived = archived[archived.method == "MLE"].copy()
        archived["method"] = "pre-repair-MLE"
        frame = pd.concat([frame, archived], ignore_index=True)
    frame["phase"] = phase
    frames.append(frame)
metrics = pd.concat(frames, ignore_index=True)
metrics["fallback_used"] = metrics.fallback.notna()
metrics["label"] = metrics.method + np.where(
    metrics.background, " + background", " / no background"
)
metrics.loc[metrics.method == "SigProfilerAssignment", "label"] = "SigProfilerAssignment"
metrics["family"] = np.select(
    [
        metrics.method == "MLE",
        metrics.method == "sparse-bic",
        metrics.method == "SigProfilerAssignment",
        metrics.method == "pre-repair-MLE",
    ],
    ["Repaired MLE", "Sparse BIC", "SigProfilerAssignment", "Pre-repair MLE"],
    default="Other configuration",
)
summary = (
    metrics.groupby(["phase", "condition", "burden", "method", "background", "label", "family"])
    .mean(numeric_only=True)
    .reset_index()
)
summary.to_csv(OUT / "summary.tsv", sep="\t", index=False)

focus = summary[
    (summary.condition == "clean")
    & summary.method.isin(["pre-repair-MLE", "MLE", "sparse-bic", "SigProfilerAssignment"])
]
for phase in ("paired", "fresh"):
    lines = [
        "| Mutations | Method | Attribution error | Active signatures | Reconstruction cosine |",
        "|---:|---|---:|---:|---:|",
    ]
    for row in focus[focus.phase == phase].sort_values(["burden", "label"]).itertuples():
        lines.append(
            f"| {row.burden} | {row.label} | {row.aggregate_null_tae:.5f} | {row.active_known:.2f} | {row.cosine:.6f} |"
        )
    (OUT / f"{phase}-clean.md").write_text("\n".join(lines) + "\n")

paired = []
rng = np.random.default_rng(202610271)
for (phase, condition, burden), group in metrics.groupby(["phase", "condition", "burden"]):
    for background in (False, True):
        for method, comparator in [
            ("MLE", "pre-repair-MLE"),
            ("sparse-bic", "MLE"),
            ("sparse-bic", "SigProfilerAssignment"),
        ]:
            first = group[(group.method == method) & (group.background == background)].set_index(
                ["cancer", "sample"]
            )
            second = group[
                (group.method == comparator)
                & (
                    group.background
                    == (False if comparator == "SigProfilerAssignment" else background)
                )
            ].set_index(["cancer", "sample"])
            difference = first.aggregate_null_tae - second.aggregate_null_tae
            assert len(difference) == 64 and not difference.isna().any()
            bootstrap = np.zeros(2000)
            for _, stratum in difference.groupby(level="cancer"):
                values = stratum.to_numpy()
                bootstrap += (
                    values[rng.integers(len(values), size=(2000, len(values)))].mean(axis=1) / 8
                )
            paired.append(
                {
                    "phase": phase,
                    "condition": condition,
                    "burden": burden,
                    "background": background,
                    "method": method,
                    "comparator": comparator,
                    "difference": difference.mean(),
                    "ci_low": np.quantile(bootstrap, 0.025),
                    "ci_high": np.quantile(bootstrap, 0.975),
                }
            )
pd.DataFrame(paired).to_csv(OUT / "paired-differences.tsv", sep="\t", index=False)

color = {
    "field": "family",
    "type": "nominal",
    "title": None,
    "scale": {
        "domain": [
            "Repaired MLE",
            "Sparse BIC",
            "SigProfilerAssignment",
            "Pre-repair MLE",
            "Other configuration",
        ],
        "range": ["#00796b", "#8e44ad", "#cf6d00", "#34495e", "#aab4be"],
    },
}
y = {
    "field": "label",
    "type": "nominal",
    "sort": "ascending",
    "title": None,
    "axis": {"labelLimit": 310},
}
chart = {
    "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
    "data": {"values": json.loads(summary.to_json(orient="records"))},
    "params": [
        {
            "name": "phase_choice",
            "value": "fresh",
            "bind": {"input": "select", "options": ["fresh", "paired"], "name": "Panel: "},
        },
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
        {
            "filter": "datum.phase === phase_choice && datum.condition === condition_choice && datum.burden === burden_choice"
        }
    ],
    "hconcat": [],
    "resolve": {"scale": {"y": "shared"}},
    "config": {
        "font": "Arial",
        "view": {"stroke": None},
        "legend": {"orient": "bottom", "columns": 3},
    },
}
for metric, title, domain in [
    ("aggregate_null_tae", "Attribution error (lower is better)", {"domain": [0, 1]}),
    ("cosine", "Reconstruction cosine (higher is better)", {"zero": False}),
]:
    chart["hconcat"].append(
        {
            "width": 290,
            "height": 600,
            "mark": {"type": "point", "filled": True, "size": 65},
            "encoding": {
                "y": y if metric == "aggregate_null_tae" else {**y, "axis": None},
                "x": {"field": metric, "type": "quantitative", "title": title, "scale": domain},
                "color": color,
                "tooltip": [
                    {"field": name, "format": ".5f"}
                    for name in [
                        "aggregate_null_tae",
                        "cosine",
                        "active_known",
                        "precision",
                        "recall",
                        "null_weight",
                        "fallback_used",
                    ]
                ],
            },
        }
    )
(FIGURES / "comparison.vl.json").write_text(json.dumps(chart, indent=2) + "\n")
(FIGURES / "comparison.svg").write_text(vlc.vegalite_to_svg(chart))
(FIGURES / "comparison.png").write_bytes(vlc.vegalite_to_png(chart, scale=1.5))
(FIGURES / "interactive.html").write_text(vlc.vegalite_to_html(chart, bundle=True))
print(focus.to_string(index=False))
