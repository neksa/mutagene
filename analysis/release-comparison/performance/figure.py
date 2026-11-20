"""Plot paired profile errors and mark the lowest mean error in each table column."""

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
LABELS = {
    "MutaGene": "MutaGene",
    "SigProfilerAssignment": "SigProfilerAssignment",
    "deconstructSigs": "deconstructSigs",
    "MuSiCal": "MuSiCal",
    "signeR": "signeR",
    "signeR+background": "signeR + background",
    "NNLS": "NNLS baseline",
}


def draw(errors, output):
    samples = pd.read_csv(HERE / "samples.tsv", sep="\t", index_col=0)
    archive = pd.read_csv(HERE.parent / "signer-panel-errors.tsv", sep="\t", index_col=0)
    if set(errors.index) != set(samples.index) or set(archive.index) != set(samples.index):
        raise ValueError("The error panels use different samples")
    errors = errors.join(archive[["signeR", "signeR+background"]], validate="one_to_one")
    errors = errors.loc[samples.index, list(LABELS)]
    if not np.isfinite(errors).all().all():
        raise ValueError("Each sample requires one finite error per method")
    means = errors.groupby(samples.burden).mean().T
    means.insert(0, "All", errors.mean())
    winners = means.eq(means.min())
    mpl.rcParams.update({"font.size": 7.6, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(7.4, 3.1))
    grid = fig.add_gridspec(
        1,
        5,
        width_ratios=[1, 1, 1, 0.08, 2.03],
        left=0.187,
        right=0.99,
        top=0.77,
        bottom=0.215,
        wspace=0.16,
    )
    axes = [fig.add_subplot(grid[0, i]) for i in (0, 1, 2, 4)]
    fig.text(0.185, 0.92, "A   Per-profile attribution error", fontsize=10)
    fig.text(axes[-1].get_position().x0, 0.92, "B   Mean attribution error", fontsize=10)
    rng = np.random.default_rng(177)
    for column, (ax, burden) in enumerate(zip(axes, (100, 2000, 50000))):
        subset = errors.loc[samples.burden == burden]
        if len(subset) != 32:
            raise ValueError(f"Expected 32 profiles with {burden} mutations")
        # Each sample keeps its vertical offset across methods within a panel.
        offsets = rng.permutation(np.linspace(-0.23, 0.23, len(subset)))
        for row, method in enumerate(LABELS):
            values = subset[method]
            color = "#056963" if method == "MutaGene" else "#52636e"
            ax.boxplot(
                values,
                positions=[row],
                orientation="horizontal",
                widths=0.50,
                manage_ticks=False,
                patch_artist=True,
                showfliers=False,
                boxprops={"facecolor": "none", "edgecolor": color, "linewidth": 0.75},
                medianprops={"color": color, "linewidth": 1.3},
                whiskerprops={"color": color, "linewidth": 0.7},
                capprops={"color": color, "linewidth": 0.7},
            )
            ax.scatter(values, row + offsets, s=3.5, color=color, alpha=0.55, linewidths=0)
        ax.set_xlim(-0.045, 1.045)
        ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
        ax.set_title(f"{burden:,} mutations", fontsize=7.6, pad=5)
        ax.set_yticks(range(len(LABELS)), list(LABELS.values()) if column == 0 else [])
        ax.grid(axis="x", color="#e1e6e5", linewidth=0.5)
        ax.set_axisbelow(True)
    table = axes[-1]
    table.set_xlim(-0.5, 3.5)
    table.set_xticks(range(4), ["All", "100", "2,000", "50,000"])
    table.xaxis.tick_top()
    table.set_yticks([])
    for row, method in enumerate(LABELS):
        for column, value in enumerate(means.loc[method]):
            marker = "*" if winners.iloc[row, column] else ""
            table.text(
                column,
                row,
                f"{value:.3f}{marker}",
                ha="center",
                va="center",
                color="#056963" if marker else "#25343b",
                fontsize=7.6,
            )
    for ax in axes:
        ax.set_ylim(len(LABELS) - 0.5, -0.5)
        ax.tick_params(axis="both", length=0, pad=3)
        ax.spines[:].set_visible(False)
    fig.text(
        0.185, 0.095, "Boxes span the middle 50%. Lines inside boxes show medians.", fontsize=7.2
    )
    fig.text(table.get_position().x0, 0.095, "* Lowest column mean.", fontsize=7.2)
    fig.text(
        0.185,
        0.025,
        "Each point represents one synthetic profile. Each burden includes 32 profiles.",
        fontsize=7.2,
    )
    fig.savefig(output, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--errors", type=Path, default=HERE / "errors.tsv")
    parser.add_argument("--output", type=Path, default=HERE / "performance.pdf")
    args = parser.parse_args()
    draw(pd.read_csv(args.errors, sep="\t", index_col=0), args.output)
