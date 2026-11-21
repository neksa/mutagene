"""Draw cohort exposures, signature spectra, motif frequencies, and driver enrichment."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import PowerNorm

HERE = Path(__file__).resolve().parent
INK = "#174c50"
TEAL = "#167e79"
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
plt.rcParams.update(
    {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42}
)


def save(fig, name, output=None):
    for target in ([output] if output else [HERE / f"{name}.pdf", HERE / f"{name}.png"]):
        fig.savefig(target, dpi=165, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


def spectra(catalog):
    profiles = pd.read_csv(HERE / "candidate-profiles.tsv", sep="\t", index_col=0)
    summary = pd.read_csv(HERE / "candidate-summary.tsv", sep="\t").set_index("experiment")
    full, control = summary.index
    near = summary.loc[full, "nearest"]
    names = ["Extra component", near, "Catalog mixture", "Recovered component", control.split()[-1]]
    panels = [
        profiles[full],
        catalog[near],
        profiles[full + " mixture"],
        profiles[control],
        catalog[names[-1]],
    ]
    inner = [f"{a}[{s}]{b}" for s in SUBS for a in "ACGT" for b in "ACGT"]
    outer = [a + b for a in "ACGT" for b in "ACGT"]
    grids = [np.array([[p[a + k + b] for k in inner] for a, b in outer]) * 100 for p in panels]
    fig = plt.figure(figsize=(8.1, 5.8))
    layout = fig.add_gridspec(
        7,
        2,
        width_ratios=[1, 0.026],
        height_ratios=[0.32, 1, 1, 1, 0.52, 1, 1],
        left=0.22,
        right=0.91,
        top=0.97,
        bottom=0.09,
        wspace=0.05,
        hspace=0.14,
    )
    for line, title in [
        (0, "A   Full catalog   ·   60 supplied signatures"),
        (4, "B   Positive control   ·   omit SBS7a_S before learning"),
    ]:
        ax = fig.add_subplot(layout[line, 0])
        ax.axis("off")
        ax.text(-0.25, 0.1, title, fontsize=13, color=INK, transform=ax.transAxes)
    norm = PowerNorm(0.5, vmin=0, vmax=max(g.max() for g in grids))
    for i, (line, name, grid) in enumerate(zip([1, 2, 3, 5, 6], names, grids)):
        ax = fig.add_subplot(layout[line, 0])
        m = ax.imshow(grid, aspect="auto", cmap="magma_r", norm=norm, interpolation="nearest")
        ax.set(yticks=[], xticks=[])
        ax.set_ylabel(name, rotation=0, ha="right", va="center", labelpad=7)
        for boundary in range(16, 96, 16):
            ax.axvline(boundary - 0.5, color="#aaa", lw=0.4)
        if i == 4:
            ax.set(
                xticks=np.arange(6) * 16 + 7.5,
                xticklabels=SUBS,
                xlabel="Inner SBS96 channels, grouped by substitution",
            )
            ax.tick_params(length=0)
    cb = fig.colorbar(m, cax=fig.add_subplot(layout[1:4, 1]))
    cb.set_label("Signature mass per channel (%)")
    save(fig, "candidate-comparison")


def motif_driver():
    x = pd.read_csv(HERE / "motifs.tsv", sep="\t")
    known = x.loc[x.family == "curated"]
    groups = ["Skin-Melanoma", "ColoRect-AdenoCA", "Lung-SCC", "Eso-AdenoCA", "Stomach-AdenoCA"]
    names = ["Melanoma", "Colorectal", "Lung SCC", "Esophageal", "Stomach"]
    fold = known.pivot(index="motif", columns="histology", values="fold")[groups]
    q = known.pivot(index="motif", columns="histology", values="q").loc[fold.index, groups]
    fig, axes = plt.subplots(
        1, 2, figsize=(9.1, 3.9), gridspec_kw={"width_ratios": [1.5, 1], "wspace": 0.49}
    )
    ax = axes[0]
    ax.imshow(np.log2(fold.values), aspect="auto", cmap="BrBG", vmin=-3, vmax=3)
    for i, j in np.ndindex(fold.shape):
        f = fold.iloc[i, j]
        ax.text(
            j,
            i,
            f"{f:.1f}" + ("*" if q.iloc[i, j] < 0.05 else ""),
            ha="center",
            va="center",
            color="white" if abs(np.log2(f)) > 1.8 else "#182d30",
            fontsize=9,
        )
    ax.set(
        xticks=range(len(groups)),
        xticklabels=names,
        yticks=range(len(fold)),
        yticklabels=fold.index,
    )
    ax.tick_params(length=0)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    ax.set_title("A   Motif enrichment (fold)", loc="left", color=INK, pad=12)
    e = pd.read_csv(HERE / "cohort-enrichment.tsv", sep="\t").set_index("unit").loc["variant"]
    ax = axes[1]
    hits = [int(e.driver_hits), int(e.other_hits)]
    totals = [int(e.driver_total), int(e.other_total)]
    rate = np.array(hits) / totals * 100
    ax.bar([0, 1], rate, width=0.57, color=[TEAL, "#acb8ba"], zorder=2)
    for j, y in enumerate(rate):
        ax.text(j, y + 0.07, f"{y:.2f}%\n{hits[j]:,}/{totals[j]:,}", ha="center", fontsize=10)
    ax.set(
        xticks=[0, 1],
        xticklabels=["Driver", "Other scored"],
        ylabel="Exact hotspot matches (%)",
        ylim=(0, max(rate) * 1.32),
    )
    ax.grid(axis="y", color="#dedede", lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title("B   BGI bladder cohort", loc="left", color=INK, pad=12)
    save(fig, "motifs-and-drivers")


def waterfall(counts, output=None):
    fits = (
        pd.read_csv(HERE / "exposures.tsv", sep="\t")
        .pivot(index="sample", columns="signature", values="exposure")
        .fillna(0)
    )
    if set(fits.index) != set(counts.columns) or not np.allclose(fits.sum(axis=1), 1, atol=0.00005):
        raise ValueError("Check sample labels and exposure totals")
    sbs3 = "SBS3_P"
    if sbs3 not in fits:
        fits[sbs3] = 0.0
    leading = list(fits.drop(columns=sbs3).mean().nlargest(5).index)
    plot = fits[[*leading, sbs3]].copy()
    plot["Other"] = fits.drop(columns=plot.columns).sum(axis=1)

    def group(sample):
        hist = sample.split("::")[0]
        if hist.startswith("Lung-"):
            return "Lung"
        return {"Skin-Melanoma": "Melanoma", "ColoRect-AdenoCA": "Colorectal"}.get(
            hist, "Other tumors"
        )

    groups = ["Melanoma", "Lung", "Colorectal", "Other tumors"]
    fig, ax = plt.subplots(figsize=(9.1, 3.45))
    colors = ["#266a76", "#e9ab5d", "#7472a7", "#cc7569", "#78a691", "#a03371", "#d9dede"]
    x = 0
    for group_name in groups:
        ordered = plot[[group(k) == group_name for k in plot.index]].sort_values(
            leading, ascending=False
        )
        pos = np.arange(len(ordered)) + x
        bottom = np.zeros(len(ordered))
        for col, color in zip(plot, colors):
            ax.bar(
                pos,
                ordered[col] * 100,
                bottom=bottom,
                width=1,
                color=color,
                linewidth=0,
                label=col if x == 0 else None,
            )
            bottom += ordered[col].to_numpy() * 100
        ax.text(pos.mean(), 104, f"{group_name} ({len(pos)})", ha="center", fontsize=10)
        x += len(ordered) + 2
    ax.set(
        ylim=(0, 100),
        xlim=(-1, x - 2),
        xticks=[],
        ylabel="Fitted exposure (%)",
        xlabel="One bar per PCAWG tumor",
    )
    ax.spines["bottom"].set_visible(False)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=7,
        frameon=False,
        fontsize=9,
        handlelength=1,
    )
    save(fig, "pcawg-exposures", output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    args = parser.parse_args()
    spectra(pd.read_csv(args.reference, sep="\t", index_col=0))
    motif_driver()
    waterfall(pd.read_csv(HERE / "pcawg-counts.tsv.gz", sep="\t", index_col=0))
