"""Measure one batch per fresh process and score each sample against the generating mixture."""

import argparse
import hashlib
import os
import subprocess
import sys
from importlib.metadata import distribution
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
METHODS = [
    "MutaGene",
    "MuSiCal",
    "NNLS",
    "MutaGene (previous)",
    "deconstructSigs",
    "SigProfilerAssignment",
]
THREADS = [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
]


def run(args):
    previous = distribution("mutagene").locate_file("mutagene/signatures/sparse.py")
    if hashlib.sha256(previous.read_bytes()).hexdigest() != (
        "1d74b3bafcd5ddb3648543e19af16ef689b77f25d5a9f1c7648daebf867bdc2c"
    ):
        raise ValueError("Install the previous MutaGene wheel specified in the benchmark README")
    folder = args.output.resolve()
    folder.mkdir(parents=True, exist_ok=False)
    samples = pd.read_csv(HERE / "samples.tsv", sep="\t", index_col=0)
    counts = pd.read_csv(HERE / "counts.tsv", sep="\t", index_col=0)
    truth = pd.read_csv(HERE / "truth.tsv", sep="\t", index_col=0)
    if args.panel == "timing":
        samples = samples.loc[samples.timing]
    counts = counts.loc[:, samples.index]
    truth = truth.loc[samples.index]
    np.testing.assert_array_equal(counts.sum(), samples.burden)
    matrix = folder / "matrix.tsv"
    counts.to_csv(matrix, sep="\t")
    reference = HERE / "reference.tsv"
    times, errors, first = [], [], {}
    rounds = 4 if args.panel == "timing" else 1
    for repeat in range(rounds):
        order = METHODS[repeat:] + METHODS[:repeat]
        for method in order:
            dest = folder / f"{repeat}-{method}"
            dest.mkdir()
            output = dest / "weights.tsv"
            env = dict(os.environ, **dict.fromkeys(THREADS, "1"))
            env.pop("PYTHONPATH", None)
            if method == "deconstructSigs":
                command = [
                    "Rscript",
                    "--vanilla",
                    HERE / "fit.R",
                    args.deconstruct.resolve(),
                    matrix,
                    reference,
                    output,
                ]
            else:
                selected = "MutaGene" if method.startswith("MutaGene") else method
                if method == "MutaGene":
                    env["PYTHONPATH"] = str(HERE.parents[2])
                command = [sys.executable, HERE / "fit.py", selected, matrix, reference, output]
            with (dest / "log.txt").open("w") as log:
                start = perf_counter()
                subprocess.run(
                    command, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, check=True
                )
                elapsed = perf_counter() - start
            fitted = pd.read_csv(output, sep="\t", index_col=0)
            if set(fitted.index) != set(truth.index) or set(fitted) != set(truth):
                raise ValueError(f"Incorrect output identifiers for {method}")
            fitted = fitted.loc[truth.index, truth.columns]
            if not np.isfinite(fitted).all().all() or fitted.min().min() < -1e-12:
                raise ValueError(f"Invalid weights for {method}")
            np.testing.assert_allclose(fitted.sum(axis=1), 1, rtol=1e-9)
            if repeat == 0:
                first[method] = fitted
                errors.append((0.5 * (truth - fitted).abs().sum(axis=1)).rename(method))
            else:
                np.testing.assert_allclose(first[method], fitted, atol=1e-12, rtol=0)
            times.append(dict(round=repeat, method=method, samples=len(samples), seconds=elapsed))
            pd.DataFrame(times).to_csv(folder / "times.tsv", sep="\t", index=False)
            print(repeat, method, round(elapsed, 3), flush=True)
    np.testing.assert_array_equal(first["MutaGene"], first["MutaGene (previous)"])
    pd.concat(errors, axis=1).to_csv(folder / "errors.tsv", sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("panel", choices=["timing", "all"])
    parser.add_argument("deconstruct", type=Path, help="Pinned deconstructSigs source directory")
    parser.add_argument("output", type=Path, help="New result directory")
    run(parser.parse_args())
