"""Plot matched candidate and reference fingerprints at two context resolutions."""

from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[1]
BASE = PAPER.parent / "analysis/release-comparison/pcawg-penta"
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
COLORS = ["#378ab5", "#41495a", "#b95e53", "#a3978c", "#58936b", "#9271a7"]


def main():
    (PAPER / "figures").mkdir(parents=True, exist_ok=True)
    reference = pd.read_csv(PAPER / "inputs/SBS1536-display.tsv", sep="\t", index_col=0)
    candidates = pd.read_csv(BASE / "candidate-profiles.tsv", sep="\t", index_col=0)
    summary = pd.read_csv(BASE / "candidate-summary.tsv", sep="\t").set_index("experiment")
    inner = [f"{a}[{sub}]{b}" for sub in SUBS for a, b in product("ACGT", repeat=2)]
    outer = ["".join(pair) for pair in product("ACGT", repeat=2)]
    fig = plt.figure(figsize=(8.0, 5.7), layout="constrained")
    subfigs = fig.subfigures(2, 1, height_ratios=[3, 2])
    for subfig, experiment, rows, title in zip(
        subfigs,
        summary.index,
        [
            [
                ("Candidate", candidates["Full catalog"]),
                ("SBS17b_P", reference.SBS17b_P),
                ("Catalog mixture", candidates["Full catalog mixture"]),
            ],
            [("Recovered", candidates["Withhold SBS7a_S"]), ("SBS7a_S", reference.SBS7a_S)],
        ],
        ["A   All 60 supplied signatures", "B   Recovery after omitting SBS7a_S"],
    ):
        axes = subfig.subplots(
            len(rows), 2, gridspec_kw={"width_ratios": [3.4, 1.25]}, squeeze=False
        )
        subfig.suptitle(title, x=0.02, ha="left", fontsize=12, color="#174c50")
        profiles = [100 * p / p.sum() for _, p in rows]
        collapsed = [p.groupby(p.index.str[1:-1]).sum().loc[inner] for p in profiles]
        center = collapsed[0].idxmax()
        details = [[p[o[0] + center + o[1]] for o in outer] for p in profiles]
        axes[0, 0].set_title("96 inner channels", fontsize=9)
        axes[0, 1].set_title(f"16 outer pairs in {center}", fontsize=9)
        for i, ((label, _), tri, detail) in enumerate(zip(rows, collapsed, details)):
            ax, zoom = axes[i]
            ax.bar(range(96), tri, color=np.repeat(COLORS, 16), width=0.87)
            zoom.bar(range(16), detail, color="#247773", width=0.82)
            ax.set_ylabel(label, rotation=0, ha="right", va="center", labelpad=8, fontsize=9)
            for a, limit in [(ax, max(s.max() for s in collapsed)), (zoom, np.max(details))]:
                a.set_ylim(0, limit * 1.08)
                a.spines[["top", "right"]].set_visible(False)
                a.tick_params(length=2, labelsize=7)
                a.set_xticks([])
            if i == len(rows) - 1:
                ax.set_xticks(np.arange(6) * 16 + 7.5, SUBS)
                zoom.set_xticks([0, 5, 10, 15], ["AA", "CC", "GG", "TT"])
    fig.savefig(PAPER / "figures/fingerprints.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
