"""Plot fitted PCAWG exposures and published signature-3 treatment subgroups."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter
from scipy.stats import binomtest

HERE = Path(__file__).resolve().parent
INK = "#174c50"
POSITIVE = "#087b76"
NEGATIVE = "#a44b2a"
DAYS_PER_MONTH = 365.25 / 12


def cohort(ax):
    weights = pd.read_csv(HERE / "pcawg-weights.tsv", sep="\t", index_col=0)
    if len(weights) != 146 or not np.allclose(weights.sum(axis=1), 1) or (weights < 0).any().any():
        raise ValueError("Expected 146 nonnegative exposure vectors with unit mass")
    leading = list(weights.drop(columns=["SBS3", "Background"]).mean().nlargest(5).index)
    grouped = weights[["SBS3", *leading]].copy()
    grouped["Other SBS"] = weights.drop(columns=["SBS3", *leading, "Background"]).sum(axis=1)
    grouped["Background"] = weights["Background"]
    histology = weights.index.str.split("::").str[0]
    groups = (
        pd.Series(histology, index=weights.index)
        .map(
            {
                "Skin-Melanoma": "Melanoma",
                "Lung-SCC": "Lung",
                "Lung-AdenoCA": "Lung",
                "ColoRect-AdenoCA": "Colorectal",
            }
        )
        .fillna("Other")
    )
    colors = [
        "#9a3c73",
        "#267b84",
        "#76bab2",
        "#ce8055",
        "#6574aa",
        "#cfaf5b",
        "#b9c1c9",
        "#454e58",
    ]
    offset = 0
    centers, labels = [], []
    for group in ["Melanoma", "Lung", "Colorectal", "Other"]:
        block = grouped.loc[groups == group].sort_values(leading, ascending=False)
        x = offset + np.arange(len(block))
        bottom = np.zeros(len(block))
        for color, name in zip(colors, grouped):
            ax.bar(
                x,
                block[name],
                width=1.01,
                bottom=bottom,
                color=color,
                linewidth=0,
                label=name if offset == 0 else None,
            )
            bottom += block[name].to_numpy()
        centers.append(x.mean())
        labels.append(f"{group}\n{len(block)}")
        offset += len(block) + 1.5
    ax.set_xlim(-0.6, x[-1] + 0.6)
    ax.set_ylim(0, 1)
    ax.set_xticks(centers, labels)
    ax.set_yticks([0, 0.5, 1])
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_ylabel("Fitted exposure")
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(-0.035, 1.04),
        ncol=4,
        frameon=False,
        columnspacing=1.5,
        handlelength=1.1,
        handletextpad=0.45,
        borderaxespad=0,
    )


def progression(ax, risk):
    records = pd.read_csv(HERE / "topacio.tsv", sep="\t")
    ticks = np.arange(0, 25, 6)
    for row, (group, color) in enumerate([("Positive", POSITIVE), ("Negative", NEGATIVE)]):
        data = records.loc[records.SigMA == group]
        table = data.pivot_table(
            index="days", columns="censor", values="count", aggfunc="sum", fill_value=0
        ).reindex(columns=[0, 1], fill_value=0)
        at_risk = table.sum(axis=1).iloc[::-1].cumsum().iloc[::-1]
        survival = (1 - table[0] / at_risk).cumprod()
        months = table.index / DAYS_PER_MONTH
        ax.step(
            [0, *months],
            [1, *survival],
            where="post",
            color=color,
            lw=1.9,
            label=f"Sig3 {group.lower()}",
        )
        censored = table[1] > 0
        ax.plot(months[censored], survival[censored], "+", color=color, ms=7, mew=1.3)
        risk.text(-1.25, row, group, ha="right", va="center", color=color, fontsize=8.5)
        for t in ticks:
            n = data.loc[data.days >= t * DAYS_PER_MONTH, "count"].sum()
            risk.text(t, row, str(n), ha="center", va="center", fontsize=8.5)
    ax.set(
        xlim=(-0.2, 24.3),
        ylim=(-0.025, 1.03),
        xticks=ticks,
        yticks=[0, 0.5, 1],
        xlabel="Months",
        ylabel="Progression-free survival",
    )
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.legend(
        frameon=False,
        ncol=2,
        handlelength=1.5,
        columnspacing=0.9,
        loc="upper right",
        bbox_to_anchor=(1.08, 1.23),
        fontsize=9,
    )
    risk.set(xlim=ax.get_xlim(), ylim=(1.55, -1))
    risk.axis("off")
    risk.text(0, -0.88, "Patients at risk", fontsize=9, ha="left")


def response(ax):
    for y, (events, size, color) in enumerate([(14, 19, POSITIVE), (4, 13, NEGATIVE)]):
        p = events / size
        low, high = binomtest(events, size).proportion_ci(method="wilson")
        ax.errorbar(p, y, xerr=[[p - low], [high - p]], fmt="o", color=color, capsize=4, ms=7)
        ax.text(p, y + 0.20, f"{events}/{size}   ({p:.0%})", ha="center", va="top", color=color)
    ax.set(
        xlim=(0, 1),
        ylim=(1.55, -0.5),
        yticks=[0, 1],
        yticklabels=["Sig3 positive", "Sig3 negative"],
        xticks=[0, 0.5, 1],
        xlabel="Objective response rate",
    )
    ax.xaxis.set_major_formatter(PercentFormatter(1))


def draw(output=None):
    plt.rcParams.update(
        {
            "font.size": 9.5,
            "pdf.fonttype": 42,
            "axes.labelcolor": "#273b41",
            "xtick.color": "#36494e",
            "ytick.color": "#36494e",
        }
    )
    fig = plt.figure(figsize=(6.8, 5.85))
    top = fig.add_axes([0.095, 0.635, 0.895, 0.205])
    km = fig.add_axes([0.11, 0.187, 0.377, 0.235])
    risk = fig.add_axes([0.11, 0.010, 0.377, 0.081])
    rate = fig.add_axes([0.665, 0.187, 0.306, 0.235])
    fig.text(0, 0.96, "A   Signature exposures in 146 PCAWG tumors", fontsize=13, color=INK)
    fig.text(0, 0.487, "B   TOPACIO · SigMA status", fontsize=13, color=INK)
    fig.text(0.568, 0.487, "C   PETREMAC · SigMA status", fontsize=13, color=INK)
    cohort(top)
    progression(km, risk)
    response(rate)
    for ax in (top, km, rate):
        ax.spines[["right", "top", "left"]].set_visible(False)
        ax.spines["bottom"].set_color("#b7c5c6")
        ax.tick_params(length=0, pad=4)
        ax.set_axisbelow(True)
    km.grid(axis="y", color="#e0e7e8", lw=0.7)
    rate.grid(axis="x", color="#e0e7e8", lw=0.7)
    for target in ([output] if output else [HERE / "cohorts.pdf", HERE / "cohorts.png"]):
        fig.savefig(target, dpi=210)
    plt.close(fig)


if __name__ == "__main__":
    draw()
