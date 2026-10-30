"""Generate paired summaries and Vega-Lite figures from the frozen benchmark outputs."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import vl_convert as vlc

ROOT = Path(__file__).resolve().parent
results = ROOT / "results"
figures = ROOT / "figures"
figures.mkdir(exist_ok=True)
m = pd.read_csv(results / "metrics.tsv", sep="\t")
labels = {
    "mle_null": "MutaGene MLE + null",
    "mle_no_null": "MutaGene MLE, no null",
    "mle_scaled_null": "Scaled MLE + null",
    "nnls": "NNLS",
    "sigprofilerassignment": "SigProfilerAssignment",
}
conditions = {
    "clean": "Clean",
    "out_of_reference_20": "20% omitted signatures",
    "uniform_20": "20% uniform noise",
    "substitution_noise_20": "20% substitution noise",
}
summary = (
    m.groupby(["condition", "burden", "method"], sort=False).mean(numeric_only=True).reset_index()
)
summary["Method"] = summary.method.map(labels)
summary["Condition"] = summary.condition.map(conditions)
summary["Mutations"] = summary.burden.astype(str)
records = json.loads(summary.to_json(orient="records"))
base = {
    "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
    "data": {"values": records},
    "facet": {"field": "Condition", "type": "nominal", "sort": list(conditions.values())},
    "columns": 2,
    "spec": {
        "width": 275,
        "height": 185,
        "mark": {"type": "line", "point": True},
        "encoding": {
            "x": {
                "field": "Mutations",
                "type": "ordinal",
                "sort": ["100", "2000", "50000"],
                "axis": {"labelAngle": 0},
                "title": "Mutations per profile",
            },
            "y": {
                "field": "aggregate_null_tae",
                "type": "quantitative",
                "title": "Aggregate-null TAE (lower is better)",
                "scale": {"domain": [0, 0.85]},
            },
            "color": {
                "field": "Method",
                "type": "nominal",
                "sort": list(labels.values()),
                "scale": {"range": ["#1565c0", "#e69f00", "#009e73", "#777777", "#cc3377"]},
            },
            "tooltip": [
                {"field": "Method"},
                {"field": "Condition"},
                {"field": "burden"},
                {"field": "aggregate_null_tae", "format": ".4f"},
                {"field": "null_weight", "format": ".4f"},
                {"field": "f1", "format": ".4f"},
            ],
        },
    },
    "config": {
        "font": "Arial",
        "axis": {"labelFontSize": 11, "titleFontSize": 11},
        "legend": {"orient": "bottom", "columns": 2, "labelLimit": 250},
        "view": {"stroke": None},
    },
}
(figures / "attribution.vl.json").write_text(json.dumps(base, indent=2) + "\n")
(figures / "attribution.svg").write_text(vlc.vegalite_to_svg(base))
(figures / "attribution.pdf").write_bytes(vlc.vegalite_to_pdf(base))
(figures / "attribution.png").write_bytes(vlc.vegalite_to_png(base, scale=2))
(figures / "interactive.html").write_text(vlc.vegalite_to_html(base, bundle=True))

rng = np.random.default_rng(20261023)
paired = []
for (condition, burden), frame in m.groupby(["condition", "burden"]):
    pivot = frame.pivot(index=["cancer", "sample"], columns="method", values="aggregate_null_tae")
    for method in labels:
        if method == "mle_null":
            continue
        difference = pivot[method] - pivot["mle_null"]
        boot = np.zeros(2000)
        for _, group in difference.groupby(level="cancer"):
            values = group.to_numpy()
            boot += values[rng.integers(0, len(values), (2000, len(values)))].mean(axis=1) / 8
        paired.append(
            {
                "condition": condition,
                "burden": int(burden),
                "method": method,
                "delta_vs_mle_null": float(difference.mean()),
                "ci_low": float(np.quantile(boot, 0.025)),
                "ci_high": float(np.quantile(boot, 0.975)),
            }
        )
pd.DataFrame(paired).to_csv(results / "paired_differences.tsv", sep="\t", index=False)
# These intervals describe Monte Carlo variation conditional on the eight sampled cancer distributions.
convergence = (
    m[m.method.str.startswith("mle")]
    .groupby(["method", "burden"])
    .optimizer_success.agg(["mean", "count"])
)
convergence.to_csv(results / "convergence.tsv", sep="\t")
summary.to_csv(results / "figure_summary.tsv", sep="\t", index=False)
print(
    summary[summary.burden == 50000][
        ["Condition", "Method", "aggregate_null_tae", "null_weight"]
    ].to_string(index=False)
)
