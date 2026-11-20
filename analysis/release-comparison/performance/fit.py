"""Fit one labelled batch and write known weights and aggregate background mass."""

import argparse
from pathlib import Path

import pandas as pd


def fit(method, matrix, reference, output):
    counts = pd.read_csv(matrix, sep="\t", index_col=0)
    catalog = pd.read_csv(reference, sep="\t", index_col=0, float_precision="round_trip")
    if not counts.index.equals(catalog.index):
        raise ValueError("Channel labels differ")
    if method == "MutaGene":
        from mutagene.signatures.identify import decompose_mutational_profile_counts

        signatures = (catalog.to_numpy(), list(catalog))
        estimates = []
        for sample in counts:
            weights, _, _ = decompose_mutational_profile_counts(
                counts[sample].to_numpy(), signatures
            )
            estimates.append(
                [*weights[: len(catalog.columns)], weights[len(catalog.columns) :].sum()]
            )
        result = pd.DataFrame(estimates, index=counts.columns, columns=[*catalog, "Background"])
    else:
        if method == "MuSiCal":
            from musical.refit import refit

            activities, _ = refit(counts, catalog)
            activities = activities.T
        elif method == "NNLS":
            from scipy.optimize import nnls

            estimates = [nnls(catalog, counts[sample])[0] for sample in counts]
            activities = pd.DataFrame(estimates, index=counts.columns, columns=catalog.columns)
        else:
            from SigProfilerAssignment import Analyzer
            from sigProfilerPlotting import get_context_reference

            order = get_context_reference("96")
            dest = output.parent
            for frame, name in ((counts, "counts"), (catalog, "reference")):
                frame.loc[order].to_csv(dest / f"{name}.tsv", sep="\t", index_label="MutationType")
            Analyzer.cosmic_fit(
                samples=str(dest / "counts.tsv"),
                output=str(dest / "assignment"),
                signature_database=str(dest / "reference.tsv"),
                input_type="matrix",
                make_plots=False,
                sample_reconstruction_plots=False,
                export_probabilities=False,
                cpu=1,
            )
            activities = pd.read_csv(
                dest
                / "assignment/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
                sep="\t",
                index_col=0,
            )
        if set(activities.index) != set(counts.columns) or not set(activities).issubset(
            catalog.columns
        ):
            raise ValueError("Unexpected sample or signature identifiers")
        result = activities.reindex(index=counts.columns, columns=catalog.columns, fill_value=0)
        result = result.div(result.sum(axis=1), axis=0)
        result["Background"] = 0.0
    result.to_csv(output, sep="\t", index_label="sample")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("method", choices=["MutaGene", "SigProfilerAssignment", "NNLS", "MuSiCal"])
    for name in ("matrix", "reference", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    fit(args.method, args.matrix, args.reference, args.output)
