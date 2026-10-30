"""Audit likelihood optimality without changing the production fitter or frozen benchmark."""

import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from run_benchmark import score, sha256
from scipy.optimize import minimize, nnls

from mutagene.signatures import get_dummy_signatures_lists
from mutagene.signatures.identify import NegLogLik, decompose_mutational_profile_counts

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "likelihood-audit"


def certificate(matrix, counts, weights):
    """Return normalized NLL and a convex simplex optimality gap upper bound."""
    weights = np.maximum(weights, 0)
    weights = weights / weights.sum()
    observed = counts > 0
    design = matrix[observed]
    target = counts[observed] / counts.sum()
    predicted = design @ weights
    if (predicted <= 0).any():
        return np.inf, np.inf
    gradient = -(design.T @ (target / predicted))
    return -target @ np.log(predicted), float(gradient @ weights - gradient.min())


def solve(matrix, counts, variant="stable"):
    """Compare likelihood parameterizations on the same full signature set."""
    observed = counts > 0
    design = matrix[observed]
    target = counts[observed] / counts.sum()
    start, _ = nnls(matrix, counts / counts.sum())
    start /= start.sum()
    start = start * (1 - 1e-6) + 1e-6 / len(start)
    factor = counts.sum() if variant == "raw-equality-gradient" else 1.0

    def objective(weights):
        predicted = np.maximum(design @ weights, 1e-300)
        return -factor * target @ np.log(predicted), -factor * design.T @ (target / predicted)

    constraint = {
        "type": "eq",
        "fun": lambda w: w.sum() - 1,
        "jac": lambda w: np.ones(len(w)),
    }
    fun, jac = objective, True
    if variant == "scaled-equality-finite-difference":
        fun, jac = lambda w: objective(w)[0], None
    elif variant == "scaled-inequality-gradient":
        constraint = {
            "type": "ineq",
            "fun": lambda w: 1 - w.sum(),
            "jac": lambda w: -np.ones(len(w)),
        }
    result = minimize(
        fun,
        start,
        jac=jac,
        method="SLSQP",
        bounds=[(0, 1)] * len(start),
        constraints=constraint,
        options={"ftol": 1e-12, "maxiter": 2000},
    )
    weights = np.maximum(result.x, 0)
    weights /= weights.sum()
    nll, gap = certificate(matrix, counts, weights)
    return weights, {
        "success": bool(result.success),
        "message": str(result.message),
        "iterations": int(result.nit),
        "nll_per_mutation": float(nll),
        "gap_bound": gap,
    }


def solve_task(task):
    return solve(*task)


def main():
    OUT.mkdir(exist_ok=True)
    source = ROOT / "sparse-results/evaluate"
    spec = json.loads((source / "specification.json").read_text())
    for filename, digest in spec["source_sha256"].items():
        assert sha256(ROOT.parents[1] / "mutagene/signatures" / filename) == digest
    with np.load(source / "inputs.npz") as archive:
        inputs = {key: archive[key] for key in archive.files}
    with np.load(source / "estimates.npz") as archive:
        estimates = {key: archive[key] for key in archive.files}
    metrics = pd.read_csv(source / "metrics.tsv", sep="\t", low_memory=False)
    metric_lookup = {
        (row.sample, row.method, row.background): row.Index for row in metrics.itertuples()
    }
    samples = pd.read_csv(source / "samples.tsv", sep="\t")
    dummy = np.array([v for _, v in get_dummy_signatures_lists()]).T
    reference = inputs["reference"]
    tasks = [
        (np.column_stack([reference, dummy]) if background else reference, counts)
        for counts in inputs["counts"]
        for background in (False, True)
    ]
    solutions = []
    with ProcessPoolExecutor(max_workers=4) as pool:
        for result in pool.map(solve_task, tasks):
            solutions.append(result)
            if len(solutions) % 128 == 0:
                print(f"Solved {len(solutions)}/1536 full MLE fits", flush=True)
    rows, fitted = [], []
    for i, sample in samples.iterrows():
        counts = inputs["counts"][i]
        for background in (False, True):
            matrix = np.column_stack([reference, dummy]) if background else reference
            weights, info = solutions[2 * i + int(background)]
            null = weights[67:] if background else np.zeros(6)
            values = score(
                weights[:67],
                null,
                inputs["truth"][i],
                counts,
                inputs["names"].tolist(),
                inputs["generating_names"].tolist(),
                reference,
                sample.null_truth,
            )
            common = {
                "sample": sample["sample"],
                "condition": sample.condition,
                "burden": sample.burden,
                "background": background,
            }
            rows.append(
                {
                    **common,
                    "method": "stable-full-MLE",
                    **info,
                    **values,
                    "active_known": int(np.count_nonzero(weights[:67] * counts.sum() >= 10)),
                }
            )
            fitted.append(np.r_[weights[:67], null])
            for method in ("MLE", "MLE-scaled", "sparse-bic", "SigProfilerAssignment"):
                if method == "SigProfilerAssignment" and background:
                    continue
                index = metric_lookup[(sample["sample"], method, background)]
                old = metrics.iloc[index]
                vector = (
                    np.r_[estimates["known"][index], estimates["null"][index]]
                    if background
                    else estimates["known"][index]
                )
                nll, gap = certificate(matrix, counts, vector)
                rows.append(
                    {
                        **common,
                        "method": method,
                        "success": old.optimizer_success,
                        "fallback": pd.notna(old.fallback),
                        "weight_sum": vector.sum(),
                        "nll_per_mutation": nll,
                        "gap_bound": gap,
                        "nll_excess_over_stable": nll - info["nll_per_mutation"],
                        "aggregate_null_tae": old.aggregate_null_tae,
                        "active_known": old.active_known,
                        "cosine": old.cosine,
                    }
                )
        if (i + 1) % 64 == 0:
            print(f"Audited {i + 1}/768 profiles", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "metrics.tsv", sep="\t", index=False)
    np.savez_compressed(OUT / "stable-estimates.npz", weights=np.array(fitted))
    frame.groupby(["condition", "burden", "background", "method"]).mean(numeric_only=True).to_csv(
        OUT / "summary.tsv", sep="\t"
    )

    # One profile per condition and burden isolates solver parameterization effects.
    ablations = []
    for _, sample in samples.groupby(["condition", "burden"]).head(1).iterrows():
        i = samples.index[samples["sample"] == sample["sample"]].item()
        for background in (False, True):
            matrix = np.column_stack([reference, dummy]) if background else reference
            for variant in (
                "stable",
                "raw-equality-gradient",
                "scaled-equality-finite-difference",
                "scaled-inequality-gradient",
            ):
                _, info = solve(matrix, inputs["counts"][i], variant)
                ablations.append(
                    {
                        "sample": sample["sample"],
                        "burden": sample.burden,
                        "background": background,
                        "variant": variant,
                        **info,
                    }
                )
    pd.DataFrame(ablations).to_csv(OUT / "ablations.tsv", sep="\t", index=False)

    # Repeating identical counts changes burden without changing the MLE optimum.
    repetition = []
    for i in (
        samples[(samples.condition == "clean") & (samples.burden == 100)]
        .groupby("cancer")
        .head(1)
        .index
    ):
        counts = inputs["counts"][i]
        matrix = np.column_stack([reference, dummy])
        stable, _ = solve(matrix, counts)
        for multiplier in (1, 20, 500):
            info = {}
            weights, _, _ = decompose_mutational_profile_counts(
                counts * multiplier,
                (reference, inputs["names"].tolist()),
                func="MLE",
                enable_dummy=True,
                diagnostics=info,
            )
            nll, gap = certificate(matrix, counts, weights)
            repetition.append(
                {
                    "sample": samples.iloc[i]["sample"],
                    "multiplier": multiplier,
                    "nll_per_mutation": nll,
                    "gap_bound": gap,
                    "exposure_l1_from_stable": np.abs(weights - stable).sum(),
                    **info,
                }
            )
    pd.DataFrame(repetition).to_csv(OUT / "count-repetition.tsv", sep="\t", index=False)

    # Closed-form binomial example checks formula, boundary behavior, and reporting units.
    matrix = np.array([[0.8, 0.2], [0.2, 0.8]])
    counts = np.array([65.0, 35.0])
    exact = np.array([0.75, 0.25])
    weights, info = solve(matrix, counts)
    assert np.allclose(weights, exact, atol=1e-6)
    epsilon = 1e-7
    right = (
        NegLogLik(exact + [epsilon, 0], matrix, counts) - NegLogLik(exact, matrix, counts)
    ) / epsilon
    left = (
        NegLogLik(exact, matrix, counts) - NegLogLik(exact - [epsilon, 0], matrix, counts)
    ) / epsilon
    reports = {}
    for method in ("MLE", "Frobenius", "Cos", "Compat"):
        fitted_weights, _, results = decompose_mutational_profile_counts(
            counts, (matrix, ["A", "B"]), func=method, enable_dummy=False
        )
        reports[method] = {
            "reported_LogLik": next(row["score"] for row in results if row["name"] == "LogLik"),
            "count_LogLik": -NegLogLik(fitted_weights, matrix, counts),
        }
    probes = {
        "binomial_exact_weights": exact.tolist(),
        "stable_weights": weights.tolist(),
        "boundary_left_derivative": left,
        "boundary_right_derivative": right,
        "impossible_channel_penalty": NegLogLik(
            np.ones(1), np.array([[1.0], [0.0]]), np.array([0, 1])
        ),
        "reporting": reports,
    }
    (OUT / "probes.json").write_text(json.dumps(probes, indent=2) + "\n")
    manifest = {
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "source_hashes": spec["source_sha256"],
        "audit_script_sha256": sha256(__file__),
        "input_hashes": {
            name: sha256(source / name)
            for name in ("inputs.npz", "estimates.npz", "metrics.tsv", "samples.tsv")
        },
        "stable_fits": len(fitted),
        "stable_successes": int(
            frame.loc[frame.method == "stable-full-MLE", "success"].astype(bool).sum()
        ),
        "maximum_stable_gap_bound": float(
            frame.loc[frame.method == "stable-full-MLE", "gap_bound"].max()
        ),
    }
    (OUT / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
