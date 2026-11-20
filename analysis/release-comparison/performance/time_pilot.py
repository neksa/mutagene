"""Repeat the original CLI timing pilot with the previous and current MutaGene fitters."""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from time import perf_counter

import pandas as pd
from benchmark import HERE, THREADS


def time_pilot(panel, output):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    matrix = panel.resolve() / "matrix.tsv"
    counts = pd.read_csv(matrix, sep="\t", index_col=0)
    if counts.shape != (96, 32):
        raise ValueError("Use the 32-profile pilot from compare.py")
    results = []
    methods = ["MutaGene", "SigProfilerAssignment", "MutaGene (previous)"]
    previous_output = None
    for repeat in range(4):
        for method in methods[repeat % 3 :] + methods[: repeat % 3]:
            dest = output / f"{repeat}-{method}"
            env = dict(os.environ, **dict.fromkeys(THREADS, "1"))
            env.pop("PYTHONPATH", None)
            if method == "SigProfilerAssignment":
                command = [
                    sys.executable,
                    HERE.parents[2] / "examples/signatures/sigprofiler_assignment.py",
                    "--matrix",
                    matrix,
                    "--reference",
                    panel.resolve() / "reference.tsv",
                    "--output",
                    dest,
                ]
            else:
                if method == "MutaGene":
                    env["PYTHONPATH"] = str(HERE.parents[2])
                command = [
                    sys.executable,
                    "-m",
                    "mutagene",
                    "signature",
                    "-f",
                    "MATRIX",
                    "-s",
                    "COSMICv3",
                    "-i",
                    matrix,
                    "-o",
                    dest,
                ]
            with dest.with_suffix(".log").open("w") as log:
                start = perf_counter()
                subprocess.run(
                    command, cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT, check=True
                )
                elapsed = perf_counter() - start
            if method == "SigProfilerAssignment":
                path = (
                    dest
                    / "assignment/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
                )
                samples = pd.read_csv(path, sep="\t", index_col=0).index
            else:
                samples = pd.read_csv(dest, sep="\t")["sample"]
                if previous_output is not None and previous_output != dest.read_bytes():
                    raise ValueError("The MutaGene versions return different exposures")
                previous_output = dest.read_bytes()
            if set(samples) != set(counts.columns):
                raise ValueError(f"Missing samples for {method}")
            results.append(
                dict(round=repeat, method=method, samples=len(counts.columns), seconds=elapsed)
            )
            pd.DataFrame(results).to_csv(output / "times.tsv", sep="\t", index=False)
            print(repeat, method, round(elapsed, 3), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("panel", type=Path, help="Pilot directory from compare.py")
    parser.add_argument("output", type=Path, help="New result directory")
    args = parser.parse_args()
    time_pilot(args.panel, args.output)
