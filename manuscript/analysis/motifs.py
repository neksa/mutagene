"""Compare PCAWG motif fractions and group bladder tumors by their APOBEC motif fraction."""

import runpy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.proportion import proportion_confint

from mutagene.io.motifs import get_known_motifs

PAPER = Path(__file__).resolve().parents[1]
SRC = PAPER.parent / "analysis/release-comparison/pcawg-penta"


def main():
    matcher = runpy.run_path(str(SRC / "motifs.py"))
    out = PAPER / "figures"
    out.mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(SRC / "pcawg-counts.tsv.gz", sep="\t", index_col=0)
    motifs = get_known_motifs()
    fractions = pd.DataFrame({m["logo"]: matcher["curated_fraction"](counts, m) for m in motifs})
    groups = fractions.index.str.split("::").str[0]
    display = groups.map(
        lambda s: {
            "Skin-Melanoma": "Melanoma",
            "ColoRect-AdenoCA": "Colorectal",
            "Lung-SCC": "Lung SCC",
            "Eso-AdenoCA": "Esophageal",
            "Stomach-AdenoCA": "Stomach",
        }.get(s, "Other")
    )
    fractions["histology"] = groups
    fractions.to_csv(out / "motif-profiles.tsv", sep="\t", index_label="tumor")
    group_names = ["Melanoma", "Lung SCC", "Colorectal", "Esophageal", "Stomach", "Other"]
    ordered = pd.concat(
        [fractions[display == s].sort_values("T[C>K]W", ascending=False) for s in group_names]
    )
    labels = [m["logo"] for m in motifs]
    fig = plt.figure(figsize=(8.3, 6.6), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[0.72, 1], width_ratios=[1.5, 1])
    ax = fig.add_subplot(gs[0, :])
    im = ax.imshow(
        ordered[labels].T.to_numpy(),
        aspect="auto",
        vmin=0,
        vmax=1,
        cmap="viridis",
        interpolation="nearest",
    )
    start = 0
    ticks = []
    names = []
    for s in group_names:
        size = (display == s).sum()
        ticks.append(start + (size - 1) / 2)
        names.append(f"{s} ({size})")
        start += size
        ax.axvline(start - 0.5, color="white", lw=0.8)
    ax.set(xticks=ticks, xticklabels=names, yticks=range(7), yticklabels=labels)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")
    ax.tick_params(length=0, labelsize=8)
    ax.set_title("A   Motif profiles of 146 PCAWG tumors", loc="left", color="#174c50")
    fig.colorbar(im, ax=ax, aspect=12, pad=0.02, label="Fraction")
    screen = matcher["screen"](counts)
    screen.to_csv(
        out / "motif-screen.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    x = screen[screen.family == "curated"]
    cols = ["Skin-Melanoma", "Lung-SCC", "ColoRect-AdenoCA", "Eso-AdenoCA", "Stomach-AdenoCA"]
    fold = x.pivot(index="motif", columns="histology", values="fold").loc[labels, cols]
    q = x.pivot(index="motif", columns="histology", values="q").loc[labels, cols]
    ax = fig.add_subplot(gs[1, 0])
    ax.imshow(
        np.log2(fold), aspect="auto", cmap="BrBG", norm=TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
    )
    for i, j in np.ndindex(fold.shape):
        v = fold.iloc[i, j]
        ax.text(
            j,
            i,
            f"{v:.1f}" + ("*" if q.iloc[i, j] < 0.05 else ""),
            ha="center",
            va="center",
            fontsize=9,
            color="white" if abs(np.log2(v)) > 1.85 else "#172d31",
        )
    ax.set(xticks=range(5), xticklabels=group_names[:5], yticks=range(7), yticklabels=labels)
    plt.setp(ax.get_xticklabels(), rotation=32, ha="right")
    ax.tick_params(length=0, labelsize=9)
    ax.set_title("B   Histology enrichment (fold)", loc="left", color="#174c50")
    bgi = pd.read_csv(PAPER / "inputs/bgi-sbs96.tsv", sep="\t", index_col=0)
    eligible = bgi[bgi.index.str.slice(2, 5).isin(["C>T", "C>G"])].sum()
    hits = bgi[bgi.index.str.fullmatch(r"T\[C>[GT]\][AT]")].sum()
    clinical = pd.read_csv(PAPER / "inputs/bgi-stage.tsv", sep="\t").set_index("SAMPLE_ID")
    assert set(eligible.index) == set(clinical.index) and clinical.index.is_unique
    clinical["eligible"] = eligible
    clinical["matches"] = hits
    clinical["fraction"] = hits / eligible.replace(0, np.nan)
    stage = clinical.TNM.str.extract(r"^T([a1-4])", expand=False)
    assert stage.notna().all()
    clinical["invasive"] = stage.isin(["2", "3", "4"])
    clinical = clinical.dropna(subset=["fraction"])
    cutoff = clinical.fraction.median()
    clinical["group"] = np.where(clinical.fraction >= cutoff, "Higher", "Lower")
    clinical.to_csv(out / "motif-stage.tsv", sep="\t")
    data = clinical.groupby("group").invasive.agg(["sum", "count"]).loc[["Lower", "Higher"]]
    table = np.column_stack([data["sum"], data["count"] - data["sum"]])
    odds, p = fisher_exact(table[::-1])
    ax = fig.add_subplot(gs[1, 1])
    for i, (group, row) in enumerate(data.iterrows()):
        lo, hi = proportion_confint(row["sum"], row["count"], method="wilson")
        point = row["sum"] / row["count"]
        ax.errorbar(
            i,
            point * 100,
            yerr=np.array([[point - lo], [hi - point]]) * 100,
            fmt="o",
            capsize=4,
            color=["#a77a56", "#0a7774"][i],
        )
        ax.text(i, hi * 100 + 3, f"{row['sum']}/{row['count']}", ha="center", fontsize=10)
    ax.set(
        xlim=(-0.55, 1.55),
        ylim=(0, 100),
        xticks=[0, 1],
        xticklabels=["Lower", "Higher"],
        ylabel="Tumors with stage T2–T4 (%)",
        xlabel="APOBEC motif fraction",
    )
    ax.set_title("C   Bladder tumor stage", loc="left", color="#174c50")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out / "motifs.pdf")
    plt.close(fig)
    print(data, cutoff, odds, p)
    # A doubled motif weight increases its conditional fraction from 1/8 to 2/9.
    rng = np.random.default_rng(908)
    rows = []
    for name, prob in [("Null", 1 / 8), ("Double motif weight", 2 / 9)]:
        controls = rng.binomial(100, 1 / 8, (200, 30))
        cases = rng.binomial(100, prob, (200, 30))
        ps = mannwhitneyu(
            cases, controls, axis=1, alternative="greater", method="asymptotic"
        ).pvalue
        rows.extend(zip([name] * len(ps), range(len(ps)), ps))
    sim = pd.DataFrame(rows, columns=["condition", "replicate", "p"])
    sim.to_csv(out / "motif-controls.tsv", sep="\t", index=False)
    print(sim.assign(detected=sim.p < 0.05).groupby("condition").detected.sum())


if __name__ == "__main__":
    main()
