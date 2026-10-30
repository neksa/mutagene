"""Fit a labelled SBS96 matrix with an explicitly selected external reference."""

import argparse
from pathlib import Path

import pandas as pd
from SigProfilerAssignment import Analyzer
from sigProfilerPlotting import get_context_reference

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--matrix", type=Path, required=True)
parser.add_argument("--reference", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
order = get_context_reference("96")
reference = pd.read_csv(args.reference, sep="\t", index_col=0)
matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
for frame in (reference, matrix):
    if not frame.index.is_unique or set(frame.index) != set(order):
        raise ValueError("Both tables must contain exactly the 96 labelled SBS channels")
    if frame.isna().any().any() or (frame < 0).any().any():
        raise ValueError("Input values must be nonnegative and present")
# SPA reorders samples internally, but its custom signature table keeps input order.
reference.loc[order].to_csv(args.output / "reference.tsv", sep="\t", index_label="MutationType")
matrix.loc[order].to_csv(args.output / "matrix.tsv", sep="\t", index_label="MutationType")
Analyzer.cosmic_fit(
    samples=str(args.output / "matrix.tsv"),
    output=str(args.output / "assignment"),
    signature_database=str(args.output / "reference.tsv"),
    input_type="matrix",
    make_plots=False,
    sample_reconstruction_plots=False,
    export_probabilities=False,
    cpu=1,
)
