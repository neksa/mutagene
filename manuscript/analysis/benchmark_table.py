"""Summarize paired attribution errors across the 12-burden benchmark."""

from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[1]
SOURCE = PAPER.parent / "analysis/release-comparison/burden-panel-means.tsv"
CONDITIONS = {
    "clean": "Clean",
    "out_of_reference_20": "Missing signatures",
    "uniform_20": "Uniform background",
    "substitution_noise_20": "Substitution-specific background",
}


def main():
    panel = pd.read_csv(SOURCE, sep="\t")
    methods = ["adaptive-BIC", "SigProfilerAssignment"]
    means = np.average(panel[methods], weights=panel.profiles, axis=0)
    declarations = [
        rf"\newcommand{{\pooledErrors}}{{{means[0]:.4f} versus {means[1]:.4f}}}",
        rf"\newcommand{{\profileCount}}{{{panel.profiles.sum():,}}}",
    ]
    rows = []
    for condition, label in CONDITIONS.items():
        group = panel.loc[panel.condition == condition]
        first, second = np.average(group[methods], weights=group.profiles, axis=0)
        rows.append(f"{label} & {first:.4f} & {second:.4f}" + r" \\")
    table = [
        r"\newcommand{\burdenTable}{\begin{tabular}{lrr}\toprule",
        r"Condition & MutaGene & SigProfilerAssignment \\ \midrule",
        *rows,
        r"\bottomrule\end{tabular}}",
    ]
    path = PAPER / "figures/burden-table.tex"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(declarations + table) + "\n")


if __name__ == "__main__":
    main()
