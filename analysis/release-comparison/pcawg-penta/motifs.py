"""Screen curated motifs and exact pentanucleotides across PCAWG histologies."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from mutagene.dna import bases_dict, complementary_nucleotide
from mutagene.io.motifs import get_known_motifs


def curated_fraction(counts, motif):
    hits, eligible = [], []
    offset = motif["position"]
    for label in counts.index:
        sequence = label[:2] + label[3] + label[-2:]
        alt = label[5]
        can_match = is_match = False
        reverse = "".join(complementary_nucleotide[base] for base in reversed(sequence))
        for seq, target in [(sequence, alt), (reverse, complementary_nucleotide[alt])]:
            if seq[2] in bases_dict[motif["ref"]] and target in bases_dict[motif["alt"]]:
                can_match = True
                segment = seq[2 - offset : 2 - offset + len(motif["motif"])]
                is_match |= all(
                    base in bases_dict[code]
                    for base, code in zip(segment, motif["motif"], strict=True)
                )
        hits.append(is_match)
        eligible.append(can_match)
    return counts.loc[hits].sum() / counts.loc[eligible].sum()


def screen(counts):
    curated = get_known_motifs()
    rows = {item["logo"]: curated_fraction(counts, item) for item in curated}
    fractions = pd.DataFrame(rows).T
    substitution = counts.index.str.slice(3, 6)
    exact = counts / counts.groupby(substitution).transform("sum")
    combined = pd.concat([fractions, exact])
    family = ["curated"] * len(fractions) + ["pentanucleotide"] * len(exact)
    reference = [item["name"].strip() for item in curated] + [""] * len(exact)
    histology = counts.columns.str.split("::").str[0]
    frames = []
    for group, size in histology.value_counts().items():
        if size < 10:
            continue
        inside = combined.loc[:, histology == group].to_numpy()
        outside = combined.loc[:, histology != group].to_numpy()
        p = mannwhitneyu(
            inside, outside, axis=1, alternative="greater", method="asymptotic", nan_policy="omit"
        ).pvalue
        frames.append(
            pd.DataFrame(
                {
                    "histology": group,
                    "n_tumors": size,
                    "n_other": outside.shape[1],
                    "family": family,
                    "motif": combined.index,
                    "annotation": reference,
                    "median_fraction": np.nanmedian(inside, axis=1),
                    "other_median_fraction": np.nanmedian(outside, axis=1),
                    "p": p,
                }
            )
        )
    result = pd.concat(frames, ignore_index=True)
    result["fold"] = result.median_fraction / result.other_median_fraction.replace(0, np.nan)
    finite = result.p.notna()
    result.loc[finite, "q"] = multipletests(result.loc[finite, "p"], method="fdr_bh")[1]
    return result.sort_values(["q", "p", "motif"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("counts", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    catalog = pd.read_csv(args.counts, sep="\t", index_col=0)
    screen(catalog).to_csv(args.output, sep="\t", index=False)
