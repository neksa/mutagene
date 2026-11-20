"""Time the two fixed-catalog signeR configurations on the published 24-profile batch."""

import argparse
import os
import subprocess
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from benchmark import HERE, THREADS

METHODS = ["signeR", "signeR+background"]


def run(output, rscript):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    samples = pd.read_csv(HERE / "samples.tsv", sep="\t", index_col=0).query("timing")
    counts = pd.read_csv(HERE / "counts.tsv", sep="\t", index_col=0)[samples.index]
    truth = pd.read_csv(HERE / "truth.tsv", sep="\t", index_col=0).loc[samples.index]
    assert counts.shape == (96, 24)
    np.testing.assert_array_equal(counts.sum(), samples.burden)
    matrix = output / "counts.tsv"
    counts.to_csv(matrix, sep="\t")
    env = dict(os.environ, **dict.fromkeys(THREADS, "1"))
    records, first = [], {}
    for repeat in range(4):
        for method in METHODS if repeat % 2 == 0 else reversed(METHODS):
            stem = output / f"{repeat}-{method}"
            weights = stem.with_suffix(".tsv")
            with stem.with_suffix(".log").open("w") as log:
                start = perf_counter()
                subprocess.run(
                    [
                        rscript,
                        "--vanilla",
                        HERE / "fit_signer.R",
                        matrix,
                        HERE / "reference.tsv",
                        HERE / "signer-background.tsv",
                        HERE / "signer-seeds.tsv",
                        method,
                        weights,
                    ],
                    cwd=output,
                    env=env,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
                seconds = perf_counter() - start
            fitted = pd.read_csv(weights, sep="\t", index_col=0)
            assert fitted.index.equals(truth.index) and fitted.columns.equals(truth.columns)
            assert np.isfinite(fitted).all().all() and fitted.min().min() >= 0
            np.testing.assert_allclose(fitted.sum(axis=1), 1, atol=1e-12, rtol=0)
            if repeat == 0:
                first[method] = fitted
            else:
                np.testing.assert_allclose(fitted, first[method], atol=1e-12, rtol=0)
            records.append(dict(round=repeat, method=method, samples=len(samples), seconds=seconds))
            pd.DataFrame(records).to_csv(output / "times.tsv", sep="\t", index=False)
            print(repeat, method, seconds, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    run(args.output, args.rscript)
