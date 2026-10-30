"""Assess the frozen validation gates and render paired attribution results."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import vl_convert as vlc

ROOT = Path(__file__).resolve().parent
METRICS = ["aggregate_null_tae", "paper_tae"]


def paired_intervals(frame, seed):
    rng = np.random.default_rng(seed)
    records = []
    for burden, group in frame.groupby("burden"):
        for metric in METRICS:
            pivot = group.pivot(
                index=["cancer", "sample_index", "condition"], columns="method", values=metric
            )
            diff = pivot["adaptive-BIC"] - pivot["joint-BIC"]
            clusters = (
                diff.groupby(["cancer", "sample_index"]).mean().unstack("sample_index").to_numpy()
            )
            assert np.isfinite(clusters).all()
            draws = np.stack(
                [rng.choice(row, (4000, len(row)), replace=True).mean(axis=1) for row in clusters]
            ).mean(axis=0)
            low, high = np.quantile(draws, [0.025, 0.975])
            records.append(
                dict(
                    burden=int(burden),
                    metric=metric,
                    difference=float(diff.mean()),
                    low=float(low),
                    high=float(high),
                )
            )
    return pd.DataFrame(records)


def main():
    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    selection = json.loads((ROOT / "selection.json").read_text())
    frame = pd.read_csv(ROOT / "validation/metrics.tsv", sep="\t")
    sensitivity = pd.read_csv(ROOT / "sensitivity/metrics.tsv", sep="\t")
    assert len(frame) == 3072 and len(sensitivity) == 512
    assert not frame.duplicated(["sample", "method"]).any()
    intervals = paired_intervals(frame, 202610294)
    intervals.to_csv(out / "paired-differences.tsv", sep="\t", index=False)
    paired_intervals(sensitivity, 202610295).to_csv(
        out / "sensitivity-differences.tsv", sep="\t", index=False
    )
    summary = frame.groupby(["burden", "method"])[METRICS + ["seconds"]].mean().reset_index()
    summary.to_csv(out / "summary.tsv", sep="\t", index=False)
    frame.groupby(["burden", "condition", "method"])[METRICS].mean().reset_index().to_csv(
        out / "conditions.tsv", sep="\t", index=False
    )
    frame.groupby(["burden", "cancer", "method"])[METRICS].mean().reset_index().to_csv(
        out / "cancers.tsv", sep="\t", index=False
    )
    sensitivity.groupby(["burden", "method"])[METRICS].mean().reset_index().to_csv(
        out / "sensitivity.tsv", sep="\t", index=False
    )
    meta = pd.read_csv(ROOT / "validation/samples.tsv", sep="\t")
    weights = np.load(ROOT / "validation/estimates.npz")["weights"].reshape(len(meta), 2, -1)
    assert np.isfinite(weights).all() and (weights >= 0).all()
    np.testing.assert_allclose(weights.sum(axis=2), 1, atol=1e-12)
    high = meta.burden.to_numpy() > selection["maximum_burden"]
    exact = np.array_equal(weights[high, 0], weights[high, 1])
    low = intervals[intervals.burden <= selection["maximum_burden"]]
    gates = dict(
        profiles=len(meta),
        primary_fits=len(frame),
        sensitivity_profiles=len(sensitivity) // 2,
        low_burden_intervals_pass=bool((low.high < 0).all()),
        high_burden_exact_match=exact,
        high_burden_profiles=int(high.sum()),
        validation_seed=202610292,
        statement="Synthetic validation with the frozen policy; patient-level generalization remains untested.",
    )
    (out / "gates.json").write_text(json.dumps(gates, indent=2) + "\n")
    text = [
        "| Mutations | Previous joint BIC | Adaptive joint BIC | Relative error reduction |",
        "|---|---:|---:|---:|",
    ]
    for burden, group in summary.groupby("burden"):
        lookup = group.set_index("method")
        old = lookup.loc["joint-BIC", "aggregate_null_tae"]
        new = lookup.loc["adaptive-BIC", "aggregate_null_tae"]
        text.append(f"| {burden:,} | {old:.5f} | {new:.5f} | {(old-new)/old:.1%} |")
    text += [
        "",
        "| Mutations | Thresholded error before | Thresholded error after |",
        "|---|---:|---:|",
    ]
    for burden, group in summary.groupby("burden"):
        lookup = group.set_index("method")
        text.append(
            f"| {burden:,} | {lookup.loc['joint-BIC','paper_tae']:.5f} | {lookup.loc['adaptive-BIC','paper_tae']:.5f} |"
        )
    (out / "tables.md").write_text("\n".join(text) + "\n")
    chart = summary.melt(
        id_vars=["burden", "method"], value_vars=METRICS, var_name="metric", value_name="error"
    )
    spec = {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "data": {"values": chart.to_dict("records")},
        "facet": {"column": {"field": "metric", "type": "nominal", "title": None}},
        "spec": {
            "width": 360,
            "height": 260,
            "mark": {"type": "line", "point": True},
            "encoding": {
                "x": {
                    "field": "burden",
                    "type": "quantitative",
                    "scale": {"type": "log", "domain": [100, 50000]},
                    "axis": {"values": [100, 250, 500, 1000, 2000, 50000]},
                    "title": "Mutation burden",
                },
                "y": {"field": "error", "type": "quantitative", "title": "Attribution error"},
                "color": {"field": "method", "type": "nominal", "title": None},
                "strokeDash": {"field": "method", "type": "nominal", "title": None},
                "tooltip": [{"field": x} for x in ["burden", "method", "error"]],
            },
        },
        "config": {"font": "Arial", "view": {"stroke": None}, "legend": {"orient": "bottom"}},
    }
    (out / "comparison.vl.json").write_text(json.dumps(spec, indent=2) + "\n")
    (out / "comparison.svg").write_text(vlc.vegalite_to_svg(spec))
    (out / "comparison.png").write_bytes(vlc.vegalite_to_png(spec, scale=2))
    (out / "interactive.html").write_text(vlc.vegalite_to_html(spec, bundle=True))
    print(json.dumps(gates, indent=2))
    print(intervals.to_string(index=False))


if __name__ == "__main__":
    main()
