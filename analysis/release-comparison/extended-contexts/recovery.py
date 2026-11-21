"""Measure exposure recovery on all pairs from six reference signatures per context."""

import argparse
import subprocess
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


def benchmark(context, workdir):
    rng = np.random.default_rng(7321)
    catalog = pd.read_csv(workdir / f"{context}-reference.tsv", sep="\t", index_col=0).iloc[:, :6]
    catalog = catalog / catalog.sum()
    reference = workdir / f"{context}-six.tsv"
    matrix = workdir / f"{context}-mixtures.tsv"
    estimates = workdir / f"{context}-estimates.tsv"
    catalog.to_csv(reference, sep="\t")
    counts, truth, records = {}, {}, []
    for burden in [100, 1000]:
        for first, second in combinations(catalog.columns, 2):
            sample = f"{first}_{second}_{burden}"
            counts[sample] = rng.multinomial(burden, 0.7 * catalog[first] + 0.3 * catalog[second])
            truth[sample] = {first: 0.7, second: 0.3}
            records.append((context, sample, burden, first, second))
    pd.DataFrame(counts, index=catalog.index).to_csv(matrix, sep="\t")
    subprocess.run(
        [
            sys.executable,
            "-m",
            "mutagene",
            "signature",
            "--context",
            context,
            "-f",
            "MATRIX",
            "--reference-matrix",
            str(reference),
            "-i",
            str(matrix),
            "-o",
            str(estimates),
        ],
        check=True,
    )
    weights = pd.read_csv(estimates, sep="\t").pivot(
        index="sample", columns="signature", values="exposure"
    )
    weights = weights.reindex(index=counts, columns=catalog.columns).fillna(0)
    if not np.allclose(weights.sum(axis=1), 1, atol=1e-5):
        raise ValueError("Every fit must have unit mass")
    target = pd.DataFrame.from_dict(truth, orient="index").fillna(0)
    errors = 0.5 * (weights - target).abs().sum(axis=1)
    results = pd.DataFrame(
        records, columns=["context", "sample", "burden", "major_signature", "minor_signature"]
    )
    results["attribution_error"] = results["sample"].map(errors)
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", type=Path)
    args = parser.parse_args()
    results = pd.concat([benchmark(context, args.workdir) for context in ("SBS1536", "ID83")])
    results.to_csv(Path(__file__).with_name("recovery.tsv"), sep="\t", index=False)
    print(results.groupby(["context", "burden"]).attribution_error.mean())
