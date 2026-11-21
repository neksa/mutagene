"""Compare exact hotspot recovery at identical review budgets."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from mutagene.io.profile import read_profile_str, read_signatures
from mutagene.mutability.enrichment import is_missense
from mutagene.mutability.mutability import (
    THRESHOLD_DRIVER,
    THRESHOLD_PASSENGER,
    calculate_base_substitution_mutability,
    calculate_codon_mutability,
    get_all_codon_substitutions,
    predict_driver,
)
from mutagene.signatures import get_dummy_signatures_lists
from mutagene.signatures.identify import decompose_mutational_profile_counts

PAPER = Path(__file__).resolve().parents[1]
DATA = PAPER.parent / "analysis/release-comparison/pcawg-penta"
FIG = PAPER / "figures"
METHODS = ["Recurrence", "Empirical B-score", "Signature B-score"]


def score_alleles():
    counts = pd.read_csv(PAPER / "inputs/bgi-sbs96.tsv", sep="\t", index_col=0)
    empirical = np.asarray(read_profile_str(counts.sum(axis=1).to_csv(sep="\t", header=False)))
    matrix, names = read_signatures("COSMICv3")
    info = {}
    weights, _, _ = decompose_mutational_profile_counts(
        empirical, (matrix, names), diagnostics=info
    )
    assert info["optimizer_success"] and info["fallback"] is None
    basis = np.column_stack(
        [matrix / matrix.sum(axis=0), *[p for _, p in get_dummy_signatures_lists()]]
    )
    reconstruction = empirical.sum() * (basis @ weights)
    assert np.isclose(reconstruction.sum(), empirical.sum())
    ranks = pd.read_csv(DATA / "bgi-ranked.tsv.gz", sep="\t", comment="#")
    alleles = ranks.loc[
        ranks.mutation.map(is_missense)
        & ranks.label.isin(["Driver", "Potential driver", "Passenger"])
    ].merge(
        pd.read_csv(PAPER / "inputs/bgi-codons.tsv.gz", sep="\t"),
        on=["gene", "mutation"],
        validate="one_to_one",
    )
    alleles["Recurrence"] = alleles.observed
    for name, profile in zip(METHODS[1:], [empirical, reconstruction]):
        model = calculate_base_substitution_mutability(profile, counts.shape[1])
        scores, probabilities = [], []
        for row in alleles.itertuples():
            substitutions = get_all_codon_substitutions(row.mutation[0], row.seq5, row.mutation[-1])
            p = calculate_codon_mutability(model, row.seq5, substitutions)
            scores.append(
                predict_driver(
                    row.observed, counts.shape[1], p, THRESHOLD_DRIVER, THRESHOLD_PASSENGER
                )[0]
            )
            probabilities.append(p)
        if name == "Empirical B-score":
            np.testing.assert_allclose(scores, alleles.bscore, rtol=5e-4, atol=0)
        alleles[name] = scores
        alleles[name + " mutability"] = probabilities
    assert len(alleles) == 6639 and alleles.in_known_variants.sum() == 87
    return alleles


def curves(table, methods, hit, cluster, budgets):
    positive = table[hit].astype(int).to_numpy()
    genes, group = np.unique(table[cluster], return_inverse=True)
    rng = np.random.default_rng(4109)
    blocks = []
    for method in methods:
        score = table[method] if method == "Recurrence" else -table[method]
        _, idx = np.unique(-score, return_inverse=True)
        blocks.append(idx)

    def count_hits(idx, weights):
        sizes = np.bincount(idx, weights)
        hits = np.bincount(idx, weights * positive)
        valid = sizes > 0
        return np.interp(budgets, np.r_[0, sizes[valid].cumsum()], np.r_[0, hits[valid].cumsum()])

    unit = np.ones(len(table))
    hits = np.array([count_hits(idx, unit) for idx in blocks])
    folds = np.divide(hits, hits[0], out=np.full_like(hits, np.nan), where=hits[0] > 0)
    rows = []
    for i, m in enumerate(methods):
        for j, k in enumerate(budgets):
            rows.append(
                (m, k, hits[i, j], hits[i, j] / positive.sum(), hits[i, j] / k, folds[i, j])
            )
    result = pd.DataFrame(
        rows,
        columns=["method", "budget", "hotspots", "recall", "precision", "relative_to_recurrence"],
    )
    # Only the tied block crossing a cutoff contributes random tie-order variation.
    for method, idx in zip(methods, blocks):
        sizes = np.bincount(idx)
        successes = np.bincount(idx, weights=positive).astype(int)
        end = sizes.cumsum()
        block = end.searchsorted(budgets)
        take = budgets - (end - sizes)[block]
        certain = (successes.cumsum() - successes)[block]
        for column, quantile in [("tie_hits_low", 0.025), ("tie_hits_high", 0.975)]:
            boundary = hypergeom.ppf(quantile, sizes[block], successes[block], take)
            result.loc[result.method == method, column] = certain + boundary
    # Resample genes together and keep all fitted scores fixed.
    intervals = []
    checkpoints = [10, 100, 1000, 2451]
    js = np.searchsorted(budgets, checkpoints)
    for _ in range(2000):
        weights = np.bincount(rng.integers(len(genes), size=len(genes)), minlength=len(genes))[
            group
        ]
        h = np.array([count_hits(idx, weights)[js] for idx in blocks])
        intervals.append((h[1:] - h[0]) / checkpoints)
    limits = np.quantile(intervals, [0.025, 0.975], axis=0)
    for i, m in enumerate(methods[1:]):
        for j, k in enumerate(checkpoints):
            sel = (result.method == m) & (result.budget == k)
            result.loc[sel, ["precision_difference_low", "precision_difference_high"]] = limits[
                :, i, j
            ]
    return result


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    alleles = score_alleles()
    known = json.loads((PAPER.parent / "mutagene/data/known_drivers.json").read_text())["genes"]
    carriers = pd.read_csv(PAPER / "inputs/bgi-carriers.tsv.gz", sep="\t").merge(
        alleles[["gene", "mutation"]], on=["gene", "mutation"], validate="many_to_one"
    )
    genes = alleles.groupby("gene").agg(
        alleles=("mutation", "count"), best=("Signature B-score", "min")
    )
    genes["Signature B-score"] = (genes.best * genes.alleles).clip(upper=1)
    genes["Recurrence"] = carriers.groupby("gene")["sample"].nunique()
    genes["known"] = genes.index.isin(known)
    genes = genes.reset_index()
    tables = []
    for level, table, methods, hit in [
        ("allele", alleles, METHODS, "in_known_variants"),
        ("gene", genes, ["Recurrence", "Signature B-score"], "known"),
    ]:
        budgets = np.unique(
            np.r_[
                np.geomspace(1, len(table), 160).astype(int),
                [10, 25, 50, 100, 250, 500, 1000, 2451, len(table)],
            ]
        )
        # Fractional recovery at tied cutoffs gives equal weight to every possible tie order.
        data = curves(table, methods, hit, "gene", budgets)
        data.insert(0, "level", level)
        tables.append(data)
    results = pd.concat(tables, ignore_index=True)
    results.to_csv(FIG / "rank-recovery.tsv", sep="\t", index=False, na_rep="NA")
    selected = results[(results.level == "allele") & results.budget.isin([10, 100, 1000, 2451])]
    table = selected.pivot(index="budget", columns="method", values="hotspots")[METHODS]
    recurrence = selected[selected.method == "Recurrence"].set_index("budget")
    tex_rows = [
        f"{k:,} & {row['Recurrence']:.1f} & "
        f"{recurrence.at[k, 'tie_hits_low']:.0f}--{recurrence.at[k, 'tie_hits_high']:.0f} & "
        f"{row['Empirical B-score']:.1f} & {row['Signature B-score']:.1f}" + r" \\"
        for k, row in table.iterrows()
    ]
    header = r"""\begin{tabular}{rrrrr}\toprule
Reviewed & Recurrence & 95\% tie range & Empirical B-score & Signature B-score \\ \midrule
"""
    (FIG / "driver-table.tex").write_text(
        header + "\n".join(tex_rows) + "\n" + r"\bottomrule\end{tabular}" + "\n"
    )
    alleles.to_csv(
        FIG / "allele-scores.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.3, 3.75), layout="constrained")
    colors = ["#bc7545", "#7376a7", "#087b76"]
    for ax, y, title in zip(
        axes,
        ["recall", "relative_to_recurrence"],
        ["A   Known hotspot recovery", "B   Enrichment over recurrence"],
    ):
        for method, color in zip(METHODS, colors):
            d = results[
                (results.level == "allele") & (results.method == method) & (results.budget >= 10)
            ]
            if y == "relative_to_recurrence" and method == "Recurrence":
                continue
            ax.plot(
                d.budget,
                d[y] * (100 if y == "recall" else 1),
                label=method,
                color=color,
                lw=2,
                ls="--" if method == "Empirical B-score" else "-",
            )
        if y == "recall":
            ax.plot(
                d.budget,
                d.budget / len(alleles) * 100,
                ":",
                color=".55",
                label="Random order (mean)",
            )
        else:
            ax.axhline(1, color=colors[0], lw=1, ls=":")
        ax.set(
            xscale="log",
            xlabel="Missense changes reviewed",
            ylabel=(
                "Observed hotspots recovered (%)"
                if y == "recall"
                else "Hotspot fraction / recurrence baseline"
            ),
        )
        ax.set_title(title, loc="left", fontsize=12, color="#174c50")
        ax.spines[["right", "top"]].set_visible(False)
        ax.grid(axis="y", alpha=0.18)
    axes[0].legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "drivers.pdf")
    plt.close(fig)
    print(results[results.budget.isin([10, 100, 1000])].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
