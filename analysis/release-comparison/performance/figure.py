"""Draw mean attribution errors and matched batch times from the recorded measurements."""

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

HERE = Path(__file__).resolve().parent
ORDER = [
    "MutaGene",
    "SigProfilerAssignment",
    "MuSiCal",
    "deconstructSigs",
    "NNLS",
    "signeR",
    "signeR+background",
]


def draw(errors, times, output):
    samples = pd.read_csv(HERE / "samples.tsv", sep="\t", index_col=0)
    archive = pd.read_csv(HERE.parent / "signer-panel-errors.tsv", sep="\t", index_col=0)
    if set(errors.index) != set(samples.index) or set(archive.index) != set(samples.index):
        raise ValueError("The error panels use different samples")
    errors = errors.join(archive[["signeR", "signeR+background"]], validate="one_to_one")[ORDER]
    groups = errors.groupby(samples.burden).mean().T
    groups.insert(0, "All", errors.mean())
    matrix = groups.to_numpy()
    colors = LinearSegmentedColormap.from_list("error", ["#f5f8f7", "#b5d2ce", "#69a4a1"])
    mpl.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(7.25, 3.0), gridspec_kw={"width_ratios": [1.12, 1]}
    )
    fig.subplots_adjust(left=0.195, right=0.975, bottom=0.23, top=0.80, wspace=0.15)
    left.set_title("A   Mean attribution error", loc="left", pad=25, fontsize=10)
    right.set_title("B   Batch wall time", loc="left", pad=25, fontsize=10)
    left.imshow(matrix, cmap=colors, vmin=0, vmax=1, aspect="auto", interpolation="nearest")
    for row, values in enumerate(matrix):
        for column, value in enumerate(values):
            left.text(column, row, f"{value:.3f}", ha="center", va="center", fontsize=7.7)
    left.set_xticks(range(4), ["All", "100", "2,000", "50,000"])
    left.xaxis.tick_top()
    left.set_yticks(
        range(len(ORDER)), [name.replace("+background", " + background") for name in ORDER]
    )
    left.set_xlabel("Mutations per profile", labelpad=6)
    for boundary in np.arange(0.5, len(ORDER), 1):
        left.axhline(boundary, color="white", linewidth=1)
    left.axvline(0.5, color="white", linewidth=2)
    right.set_xscale("log")
    right.set_xlim(0.35, 160)
    right.set_xticks([0.5, 1, 3, 10, 30, 100], ["0.5", "1", "3", "10", "30", "100"])
    right.set_xlabel("Seconds per 24 profiles (log scale)")
    right.set_ylim(left.get_ylim())
    right.set_yticks([])
    right.grid(axis="x", color="#e3e8e7", linewidth=0.6)
    right.set_axisbelow(True)
    for row, method in enumerate(ORDER):
        values = times.loc[(times.method == method) & (times["round"] > 0), "seconds"]
        if method.startswith("signeR"):
            right.text(0.43, row, "No matched timing", color="#657270", va="center", fontsize=7)
            continue
        if len(values) != 3:
            raise ValueError(f"Expected three measured times for {method}")
        median = values.median()
        color = "#106661" if method == "MutaGene" else "#345057"
        right.plot([values.min(), values.max()], [row, row], color=color, linewidth=1.4)
        right.scatter(median, row, color=color, s=22, zorder=3)
        right.annotate(
            f"{median:.2f}",
            (values.max(), row),
            xytext=(5, 0),
            textcoords="offset points",
            va="center",
            fontsize=7.5,
            color=color,
        )
    for axis in (left, right):
        axis.spines[:].set_visible(False)
        axis.tick_params(which="both", length=0)
    fig.text(
        0.018,
        0.027,
        "Error uses 96 synthetic profiles, 32 per burden. Times show the median and range of three runs on one Apple M1 Max.",
        fontsize=7.2,
        color="#405452",
    )
    fig.savefig(output, metadata={"Creator": "MutaGene benchmark"})
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in (
        ("errors", "errors.tsv"),
        ("times", "times.tsv"),
        ("output", "performance.pdf"),
    ):
        parser.add_argument(f"--{name}", type=Path, default=HERE / default)
    args = parser.parse_args()
    draw(
        pd.read_csv(args.errors, sep="\t", index_col=0),
        pd.read_csv(args.times, sep="\t"),
        args.output,
    )
