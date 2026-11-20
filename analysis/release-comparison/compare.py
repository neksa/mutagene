"""Compare attribution on paired clean mixtures from the bundled COSMICv3 catalog."""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from mutagene.io.matrix import channel_labels
from mutagene.io.profile import read_signatures
from mutagene.signatures.identify import decompose_mutational_profile_counts


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    matrix, names = read_signatures("COSMICv3")
    matrix = matrix / matrix.sum(axis=0)
    catalog = [name.removeprefix("COSMICv3-") for name in names]
    reference = pd.DataFrame(matrix, index=channel_labels(), columns=catalog)
    reference.to_csv(output / "reference.tsv", sep="\t", index_label="Type")

    rng = np.random.default_rng(20261120)
    counts = {}
    truths = {}
    mixtures = {}
    for pair in range(16):
        first, second = rng.choice(len(names), 2, replace=False)
        truth = np.zeros(len(names))
        truth[[first, second]] = [0.7, 0.3]
        for burden in (100, 2000):
            sample = f"pair{pair + 1:02d}_{burden}"
            counts[sample] = rng.multinomial(burden, matrix @ truth)
            truths[sample] = truth
            mixtures[sample] = [catalog[first], catalog[second]]

    pd.DataFrame(counts, index=reference.index).to_csv(
        output / "matrix.tsv", sep="\t", index_label="Type"
    )
    adapter = Path(__file__).resolve().parents[2] / "examples/signatures/sigprofiler_assignment.py"
    subprocess.run(
        [
            sys.executable,
            str(adapter),
            "--matrix",
            str(output / "matrix.tsv"),
            "--reference",
            str(output / "reference.tsv"),
            "--output",
            str(output / "spa"),
        ],
        check=True,
    )
    activities = pd.read_csv(
        output / "spa/assignment/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
        sep="\t",
        index_col=0,
    )
    if set(activities.index) != set(counts) or set(activities.columns) != set(catalog):
        raise ValueError("SigProfilerAssignment returns unexpected sample or signature labels")

    rows = []
    for sample, observed in counts.items():
        burden = int(observed.sum())
        fitted, _, _ = decompose_mutational_profile_counts(observed, (matrix, names))
        truth = truths[sample]
        spa = activities.loc[sample, catalog].to_numpy() / burden
        np.testing.assert_allclose(spa.sum(), 1, atol=1 / burden)
        mg_error = (abs(fitted[: len(names)] - truth).sum() + fitted[len(names) :].sum()) / 2
        spa_error = abs(spa - truth).sum() / 2
        rows.append([sample, burden, *mixtures[sample], mg_error, spa_error])
    scores = pd.DataFrame(
        rows, columns=["sample", "burden", "signature_70", "signature_30", "MutaGene", "SPA"]
    )
    scores.to_csv(output / "scores.tsv", sep="\t", index=False)
    print(scores.groupby("burden")[["MutaGene", "SPA"]].mean())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New directory for inputs, estimates, and scores")
    main(parser.parse_args().output.resolve())
