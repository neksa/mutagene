"""Render matched signeR summaries with explicit sampling limitations."""

import json
from pathlib import Path

import pandas as pd
import vl_convert as vlc

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "signer-results"
summary = pd.read_csv(OUT / "summary.tsv", sep="\t")
overall = pd.read_csv(OUT / "overall.tsv", sep="\t")
text = [
    "### All conditions and burdens",
    "",
    "| Method | Total attribution error | Thresholded known error | Active known | Background error |",
    "|---|---:|---:|---:|---:|",
]
for row in overall.sort_values("aggregate_null_tae").itertuples():
    text.append(
        f"| {row.method} | {row.aggregate_null_tae:.5f} | {row.paper_tae:.5f} | {row.active_known:.2f} | {row.null_absolute_error:.5f} |"
    )
for burden, group in summary[summary.condition == "clean"].groupby("burden"):
    text += [
        "",
        f"### Clean profiles: {burden:,} mutations",
        "",
        "| Method | Total attribution error | Background mass | Active known | Cosine |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in group.sort_values("aggregate_null_tae").itertuples():
        text.append(
            f"| {row.method} | {row.aggregate_null_tae:.5f} | {row.null_weight:.5f} | {row.active_known:.2f} | {row.cosine:.5f} |"
        )
(OUT / "tables.md").write_text("\n".join(text) + "\n")
spec = {
    "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
    "data": {"values": summary.to_dict(orient="records")},
    "facet": {
        "row": {"field": "condition", "type": "nominal", "title": None},
        "column": {
            "field": "burden",
            "type": "ordinal",
            "sort": [100, 2000, 50000],
            "title": "Mutation burden",
        },
    },
    "spec": {
        "width": 240,
        "height": 190,
        "mark": "bar",
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
                {"field": name}
                for name in ["method", "aggregate_null_tae", "paper_tae", "null_weight", "cosine"]
            ],
        },
    },
    "resolve": {"scale": {"x": "independent"}},
    "config": {"font": "Arial", "view": {"stroke": None}},
}
(OUT / "comparison.vl.json").write_text(json.dumps(spec, indent=2) + "\n")
(OUT / "comparison.svg").write_text(vlc.vegalite_to_svg(spec))
(OUT / "comparison.png").write_bytes(vlc.vegalite_to_png(spec, scale=2))
(OUT / "interactive.html").write_text(vlc.vegalite_to_html(spec, bundle=True))
print(overall.to_string(index=False))
