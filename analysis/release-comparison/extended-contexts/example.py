"""Fit published SBS1536 and ID83 examples through the MutaGene CLI."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = "https://raw.githubusercontent.com/"
HERE = Path(__file__).resolve().parent
SOURCES = {
    "SBS1536": (
        "getzlab/SignatureAnalyzer/f205be05f6e180bd04c3beaffc354b73ae3ccf9f/signatureanalyzer/ref/PCAWG/sa_PCAWG_sbs.tsv",
        "SigProfilerSuite/SigProfilerMatrixGenerator/7bbd3dd5dfff7d90de28e091cc2e2c98472f1512/SigProfilerMatrixGenerator/references/matrix/BRCA_bench/SBS/BRCA_bench.SBS1536.all",
    ),
    "ID83": (
        "SigProfilerSuite/SigProfilerAssignment/308ea01dce1ba1136aac6f2246f535284b5234d8/SigProfilerAssignment/data/Reference_Signatures/GRCh37/COSMIC_v3.4_ID_GRCh37.txt",
        "SigProfilerSuite/Protocol_SigProfiler/25e7b10bb65b871363a25ebda69fe7101c91e578/input_example/matrices/LUAD.ID83.all",
    ),
}


def run(workdir, output):
    workdir.mkdir(parents=True, exist_ok=True)
    fits, leading = [], {}
    for context, (reference_url, counts_url) in SOURCES.items():
        if context == "SBS1536":
            reference = pd.read_csv(ROOT + reference_url, sep="\t").set_index(
                "Somatic Mutation Type"
            )
            reference = reference.drop(columns=["feature", "Substitution Type"])
        else:
            reference = pd.read_csv(ROOT + reference_url, sep="\t", index_col=0)
        counts = pd.read_csv(ROOT + counts_url, sep="\t", index_col=0).iloc[:, :1]
        catalog = workdir / f"{context}-reference.tsv"
        matrix = workdir / f"{context}-counts.tsv"
        reference.to_csv(catalog, sep="\t", index_label="Type")
        counts.to_csv(matrix, sep="\t", index_label="MutationType")
        fit = workdir / f"{context}-fit.tsv"
        diagnostics = workdir / f"{context}-diagnostics.json"
        subprocess.run(
            [
                sys.executable,
                "-m",
                "mutagene",
                "signature",
                "-f",
                "MATRIX",
                "--context",
                context,
                "--reference-matrix",
                str(catalog),
                "-i",
                str(matrix),
                "-o",
                str(fit),
                "--diagnostics-out",
                str(diagnostics),
            ],
            check=True,
        )
        result = pd.read_csv(fit, sep="\t")
        info = json.loads(diagnostics.read_text())["samples"][counts.columns[0]]
        if info["fallback"] is not None or info["mutations"] != counts.sum().iloc[0]:
            raise ValueError("The fit requires review")
        result.insert(0, "context", context)
        result["sample_mutations"] = info["mutations"]
        fits.append(result)
        lead = result.loc[result.exposure.idxmax(), "signature"]
        leading[context] = reference[lead] / reference[lead].sum()
    pd.concat(fits).to_csv(output / "example-fits.tsv", index=False, sep="\t")
    draw(leading["SBS1536"], leading["ID83"], output)


def draw(penta, indel, output):
    # Rows specify the two outer bases. Columns specify the inner SBS96 channel.
    bases = "ACGT"
    substitutions = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    inner = [f"{a}[{sub}]{b}" for sub in substitutions for a in bases for b in bases]
    flanks = [a + b for a in bases for b in bases]
    grid = np.array([[penta[a + center + b] for center in inner] for a, b in flanks]) * 100
    col = grid.sum(axis=0).argmax()
    detail = grid[:, col].reshape(4, 4)
    plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42})
    fig = plt.figure(figsize=(6.8, 4.15))
    fig.text(
        0.01, 0.95, f"A   Pentanucleotide signature {penta.name}", color="#174c50", fontsize=13
    )
    heat = fig.add_axes([0.10, 0.58, 0.73, 0.28])
    colorbar = fig.add_axes([0.856, 0.58, 0.017, 0.28])
    p = heat.imshow(grid, cmap="magma_r", aspect="auto", interpolation="nearest", vmin=0)
    fig.colorbar(p, cax=colorbar, label="Mass (%)")
    heat.set(
        xticks=np.arange(6) * 16 + 7.5,
        xticklabels=substitutions,
        yticks=[0, 5, 10, 15],
        yticklabels=[flanks[i] for i in [0, 5, 10, 15]],
        ylabel="Outer bases",
    )
    for edge in np.arange(1, 6) * 16 - 0.5:
        heat.axvline(edge, color="#ccc", lw=0.45)
    fig.text(0.01, 0.44, f"B   Within {inner[col]}", color="#174c50", fontsize=12)
    ax = fig.add_axes([0.11, 0.12, 0.245, 0.25])
    ax.imshow(detail, cmap="magma_r", vmin=0, vmax=grid.max(), aspect="auto")
    ax.set(
        xticks=range(4),
        xticklabels=list(bases),
        yticks=range(4),
        yticklabels=list(bases),
        xlabel="3′ outer base",
        ylabel="5′ outer base",
    )
    for i in range(4):
        for j in range(4):
            ax.text(
                j,
                i,
                f"{detail[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=7,
                color="white" if detail[i, j] > 0.55 * grid.max() else "#21262a",
            )
    fig.text(0.46, 0.44, f"C   Indel signature {indel.name}", color="#174c50", fontsize=12)
    ax = fig.add_axes([0.62, 0.12, 0.35, 0.25])
    top = indel.nlargest(5).sort_values()
    ax.barh(top.index, top * 100, color="#39858a", height=0.62)
    ax.set_xlabel("Signature mass (%)")
    ax.spines[:].set_visible(False)
    ax.tick_params(length=0)
    for ext in ["png", "pdf"]:
        fig.savefig(output / f"contexts.{ext}", dpi=210)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("--output", type=Path, default=HERE)
    args = parser.parse_args()
    run(args.workdir, args.output)
