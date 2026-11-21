"""Compare the five PCAWG components with independently extracted SBS96 profiles."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from sklearn.metrics.pairwise import cosine_similarity

HERE = Path(__file__).resolve().parent
COMMIT = "673737e6aa492c6f31878c6486ceed13f338466a"
SOURCE = f"https://raw.githubusercontent.com/Nik-Zainal-Group/signature.tools.lib/{COMMIT}/data/RefSigSBS_v2.03"


def compare(bgi, output):
    candidates = pd.read_csv(HERE / "extracted-signatures.tsv", sep="\t", index_col=0)
    profiles = {"BGI": pd.read_csv(bgi, sep="\t", index_col=0)}
    for cohort in ["GEL", "Hartwig"]:
        url = f"{SOURCE}/OrganSpecificSigs_{cohort}_SBS_v2.03.tsv"
        profiles[cohort] = pd.read_csv(url, sep="\t", index_col=0)
    profiles["COSMIC v3"] = pd.read_csv(
        HERE.parent / "performance/reference.tsv", sep="\t", index_col=0
    )
    rows = []
    focus = {}
    for cohort, catalog in profiles.items():
        if set(catalog.index) != set(candidates.index) or not catalog.index.is_unique:
            raise ValueError(f"Unexpected SBS96 channels in {cohort}")
        catalog = catalog.loc[candidates.index]
        scores = cosine_similarity(candidates.T, catalog.T)
        for name, values in zip(candidates, scores):
            nearest = catalog.columns[values.argmax()]
            rows.append((name, cohort, catalog.shape[1], nearest, values.max()))
            if name == "Component 3":
                focus[cohort] = catalog[nearest] / catalog[nearest].sum()
    result = pd.DataFrame(
        rows, columns=["component", "cohort", "profiles_compared", "nearest", "cosine"]
    )
    result.to_csv(output / "replication.tsv", sep="\t", index=False)
    basis = profiles["COSMIC v3"].loc[candidates.index, ["SBS2", "SBS13"]]
    basis = basis / basis.sum()
    weights, _ = nnls(basis, candidates["Component 3"])
    weights /= weights.sum()
    mixture = basis @ weights
    print(
        "Component 3, SBS2/SBS13 weights and cosine:",
        weights,
        cosine_similarity([mixture], [candidates["Component 3"]])[0, 0],
    )
    draw(candidates, result, focus, output)


def draw(candidates, result, focus, output):
    colors = ["#4397b2", "#374747", "#cd755e", "#b8a28b", "#78aa82", "#9881b2"]
    types = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    order = [f"{left}[{sub}]{right}" for sub in types for left in "ACGT" for right in "ACGT"]
    plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42})
    fig, axes = plt.subplots(4, 1, figsize=(6.8, 3.6), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.235, right=0.99, top=0.86, bottom=0.12, hspace=0.20)
    fig.text(
        0,
        0.955,
        "PCAWG component 3 and its closest external profiles",
        fontsize=13,
        color="#174c50",
    )
    entries = [
        ("PCAWG", candidates["Component 3"]),
        ("BGI", focus["BGI"]),
        ("GEL", focus["GEL"]),
        ("Hartwig", focus["Hartwig"]),
    ]
    for ax, (cohort, row) in zip(axes, entries):
        ax.bar(np.arange(96), 100 * row.loc[order], color=np.repeat(colors, 16), width=0.90)
        ax.spines[:].set_visible(False)
        ax.set(
            xlim=(-1, 96), yticks=[0, 20], ylim=(0, max(25, 100 * max(r.max() for _, r in entries)))
        )
        ax.tick_params(length=0, labelsize=7)
        ax.axhline(0, color="#aab8ba", lw=0.6)
        name = "Component 3"
        if cohort != "PCAWG":
            hit = result[(result.component == name) & (result.cohort == cohort)].iloc[0]
            name = hit.nearest.split("_")[-1] if cohort != "BGI" else hit.nearest
            name += f" · {hit.cosine:.3f}"
        ax.text(-0.065, 0.7, cohort, transform=ax.transAxes, ha="right", color="#174c50")
        ax.text(-0.065, 0.28, name, transform=ax.transAxes, ha="right", fontsize=8)
    axes[-1].set_xticks(7.5 + 16 * np.arange(6), types)
    fig.text(
        0.20,
        0.011,
        "Height gives percent of signature mass per SBS96 channel. Labels give cosine similarity.",
        fontsize=8,
    )
    for ext in ["pdf", "png"]:
        fig.savefig(output / f"replication.{ext}", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bgi", type=Path)
    parser.add_argument("--output", type=Path, default=HERE)
    args = parser.parse_args()
    compare(args.bgi, args.output)
