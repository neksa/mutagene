"""Summarize complete background-policy panels without changing selection rules."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import vl_convert as vlc

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "background-results"
frames = []
intervals = []
rng = np.random.default_rng(202610275)
for phase in ("paired", "fresh"):
    frame = pd.read_csv(OUT / phase / "metrics.tsv", sep="\t", low_memory=False)
    assert len(frame) == 4608 and frame.error.isna().all()
    frame["phase"] = phase
    frames.append(frame)
    for (condition, burden), group in frame.groupby(["condition", "burden"]):
        for metric in ("aggregate_null_tae", "paper_tae"):
            pivot = group.pivot(index=["sample", "cancer"], columns="method", values=metric)
            for candidate in ("joint-BIC", "block-BIC"):
                for baseline in ("MLE", "fixed-BIC", "no-background-BIC", "SigProfilerAssignment"):
                    difference = pivot[candidate] - pivot[baseline]
                    strata = [chunk.to_numpy() for _, chunk in difference.groupby(level="cancer")]
                    draws = np.concatenate(
                        [rng.choice(chunk, size=(2000, len(chunk))) for chunk in strata], axis=1
                    ).mean(axis=1)
                    low, high = np.quantile(draws, [0.025, 0.975])
                    intervals.append(
                        dict(
                            metric=metric,
                            phase=phase,
                            condition=condition,
                            burden=burden,
                            candidate=candidate,
                            baseline=baseline,
                            difference=difference.mean(),
                            low=low,
                            high=high,
                        )
                    )
frame = pd.concat(frames, ignore_index=True)
summary = (
    frame.groupby(["phase", "condition", "burden", "method"])[
        [
            "paper_tae",
            "aggregate_null_tae",
            "active_known",
            "null_weight",
            "null_absolute_error",
            "cosine",
            "seconds",
        ]
    ]
    .mean()
    .reset_index()
)
summary.to_csv(OUT / "summary.tsv", sep="\t", index=False)
pd.DataFrame(intervals).to_csv(OUT / "paired-differences.tsv", sep="\t", index=False)
for phase in ("paired", "fresh"):
    selected = summary[(summary.phase == phase) & (summary.condition == "clean")]
    text = []
    for burden, group in selected.groupby("burden"):
        text += [
            f"### {burden:,} mutations",
            "",
            "| Method | Total attribution error | Thresholded known error | Background mass | Active known | Cosine |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        for row in group.sort_values("aggregate_null_tae").itertuples():
            text.append(
                f"| {row.method} | {row.aggregate_null_tae:.5f} | {row.paper_tae:.5f} | {row.null_weight:.5f} | {row.active_known:.2f} | {row.cosine:.5f} |"
            )
        text.append("")
    (OUT / f"{phase}-clean.md").write_text("\n".join(text))
chart_data = summary[summary.phase == "fresh"].to_dict(orient="records")
spec = {
    "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
    "description": "Fresh simulated profiles: background selection attribution error by burden and condition.",
    "data": {"values": chart_data},
    "facet": {
        "row": {"field": "condition", "type": "nominal", "title": None},
        "column": {
            "field": "burden",
            "type": "ordinal",
            "title": "Mutation burden",
            "sort": [100, 2000, 50000],
        },
    },
    "spec": {
        "width": 240,
        "height": 150,
        "mark": {"type": "bar"},
        "encoding": {
            "y": {"field": "method", "type": "nominal", "title": None},
            "x": {
                "field": "aggregate_null_tae",
                "type": "quantitative",
                "title": "Total attribution error",
                "axis": {"format": ".2~f"},
            },
            "color": {"field": "method", "type": "nominal", "legend": None},
            "tooltip": [
                {"field": x}
                for x in ["method", "aggregate_null_tae", "paper_tae", "null_weight", "cosine"]
            ],
        },
    },
    "resolve": {"scale": {"x": "independent"}},
    "config": {"font": "Arial", "axis": {"labelFontSize": 11}, "view": {"stroke": None}},
}
(OUT / "comparison.vl.json").write_text(json.dumps(spec, indent=2) + "\n")
(OUT / "comparison.svg").write_text(vlc.vegalite_to_svg(spec))
(OUT / "comparison.png").write_bytes(vlc.vegalite_to_png(spec, scale=2))
(OUT / "interactive.html").write_text(vlc.vegalite_to_html(spec, bundle=True))
print(summary[(summary.phase == "fresh") & (summary.condition == "clean")].to_string(index=False))
