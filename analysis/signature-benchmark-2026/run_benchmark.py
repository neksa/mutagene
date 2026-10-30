"""Assess MutaGene with the SigFitTest design and explicit null-signature controls."""

import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing
import platform
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import frozen_identify as identify  # noqa: E402

from mutagene.io.profile import get_profile_attributes_dict  # noqa: E402
from mutagene.signatures import get_dummy_signatures_lists  # noqa: E402

CANCERS = (
    "Liver-HCC",
    "Stomach-AdenoCA",
    "Head-SCC",
    "ColoRect-AdenoCA",
    "Lung-AdenoCA",
    "Skin-Melanoma",
    "Lymph-BNHL",
    "CNS-GBM",
)
OUT_SIGNATURES = (
    "SBS10c",
    "SBS10d",
    "SBS86",
    "SBS87",
    "SBS88",
    "SBS89",
    "SBS90",
    "SBS91",
    "SBS92",
    "SBS93",
    "SBS94",
    "SBS95",
)
CONDITIONS = ("clean", "out_of_reference_20", "uniform_20", "substitution_noise_20")
METHODS = ("mle_null", "mle_no_null", "nnls", "mle_scaled_null")


def channel_order():
    return [
        f"{a['context'][0]}[{a['mutation'][0]}>{a['mutation'][1]}]{a['context'][1]}"
        for a in get_profile_attributes_dict()
    ]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def generate(source, samples, replicates, seed, cancers):
    """Use empirical weights, the ten-mutation truth cutoff, and categorical sampling."""
    order = channel_order()
    reference = pd.read_csv(source / "input/COSMIC_v3_SBS_GRCh38.txt", sep="\t", index_col=0).loc[
        order
    ]
    generating = pd.read_csv(
        source / "input/COSMIC_v3.3.1_SBS_GRCh38.txt", sep="\t", index_col=0
    ).loc[order]
    for matrix in (reference, generating):
        if matrix.shape[0] != 96 or not matrix.index.is_unique:
            raise ValueError("Reference channels must contain 96 unique labels")
        matrix /= matrix.sum(axis=0)
    dummy = np.array([v for _, v in get_dummy_signatures_lists()]).T
    metadata, catalogs, truths = [], [], []
    for cancer_index, cancer in enumerate(cancers):
        empirical = pd.read_csv(
            source / f"cosmic tissue data/signature_contributions_WGS_{cancer}.dat",
            sep="\t",
            index_col=0,
        )
        for replicate in range(replicates):
            rng = np.random.default_rng(np.random.SeedSequence([seed, cancer_index, replicate]))
            selected = empirical.iloc[rng.choice(len(empirical), samples, replace=True)].to_numpy(
                dtype=float
            )
            for burden in (100, 2000, 50000):
                weights = selected / selected.sum(axis=1, keepdims=True) * burden
                weights[weights < 10] = 0
                if np.any(weights.sum(axis=1) == 0):
                    raise ValueError("Truth cutoff removes every signature from a sample")
                weights /= weights.sum(axis=1, keepdims=True)
                for condition in CONDITIONS:
                    for sample_index, row in enumerate(weights):
                        truth = pd.Series(0.0, index=generating.columns)
                        truth.loc[empirical.columns] = row
                        null_truth = 0.0
                        oor = []
                        if condition != "clean":
                            truth *= 0.8
                            null_truth = 0.2
                        if condition == "out_of_reference_20":
                            oor = rng.choice(OUT_SIGNATURES, 2, replace=False).tolist()
                            truth.loc[oor] = 0.1
                        probabilities = generating.to_numpy() @ truth.to_numpy()
                        if condition == "uniform_20":
                            probabilities += 0.2 / 96
                        elif condition == "substitution_noise_20":
                            probabilities += (
                                0.2 * dummy @ np.array([0.05, 0.1, 0.6, 0.1, 0.1, 0.05])
                            )
                        probabilities /= probabilities.sum()
                        counts = np.bincount(
                            rng.choice(96, size=burden, p=probabilities), minlength=96
                        )
                        identifier = f"{cancer}_{replicate}_{burden}_{condition}_{sample_index}"
                        metadata.append(
                            {
                                "sample": identifier,
                                "cancer": cancer,
                                "replicate": replicate,
                                "burden": burden,
                                "condition": condition,
                                "sample_index": sample_index,
                                "null_truth": null_truth,
                                "out_signatures": ",".join(oor),
                            }
                        )
                        catalogs.append(counts)
                        truths.append(truth.to_numpy())
    return reference, generating, np.array(catalogs), np.array(truths), pd.DataFrame(metadata)


def fit_one(task):
    index, counts, matrix, names = task
    output = []
    original_minimize = identify.minimize
    for method in METHODS:
        diagnostics = {}
        start = time.perf_counter()
        if method == "nnls":
            known, _ = nnls(matrix, counts.astype(float))
            known /= counts.sum()
            null = np.zeros(6)
        else:

            def capture(fun, x0, *args, **kwargs):
                if method == "mle_scaled_null":

                    def objective(x, a, b):
                        return fun(x, a, b) / float(np.sum(b))

                else:
                    objective = fun
                result = original_minimize(objective, x0, *args, **kwargs)
                diagnostics.update(
                    {
                        "optimizer_success": bool(result.success),
                        "optimizer_message": str(result.message),
                        "iterations": int(result.nit),
                    }
                )
                return result

            identify.minimize = capture
            try:
                h, _, _ = identify.decompose_mutational_profile_counts(
                    counts,
                    (matrix, names),
                    func="MLE",
                    others_threshold=0.0,
                    enable_dummy=method != "mle_no_null",
                )
            finally:
                identify.minimize = original_minimize
            known = h[: len(names)]
            null = h[len(names) :] if method != "mle_no_null" else np.zeros(6)
        output.append((index, method, known, null, time.perf_counter() - start, diagnostics))
    return output


def score(known, null, truth, counts, names, generating_names, matrix, null_truth):
    """Preserve unassigned mass and report both paper and aggregate-null errors."""
    burden = counts.sum()
    positions = [generating_names.index(n) for n in names]
    target = truth[positions]
    thresholded = known.copy()
    thresholded[thresholded * burden < 10] = 0
    expanded = np.zeros(len(truth))
    expanded[positions] = thresholded
    tp = int(np.count_nonzero((thresholded > 0) & (target > 0)))
    fp = int(np.count_nonzero((thresholded > 0) & (target == 0)))
    fn = int(np.count_nonzero((thresholded == 0) & (target > 0)))
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    dummy = np.array([v for _, v in get_dummy_signatures_lists()]).T
    reconstructed = matrix @ known + dummy @ null
    observed = counts / burden
    norm = np.linalg.norm(reconstructed) * np.linalg.norm(observed)
    return {
        "paper_tae": float(np.abs(expanded - truth).sum() / 2),
        "known_l1": float(np.abs(known - target).sum()),
        "aggregate_null_tae": float(
            (np.abs(known - target).sum() + abs(null.sum() - null_truth)) / 2
        ),
        "false_positive_weight": float(thresholded[target == 0].sum()),
        "false_negative_weight": float(target[thresholded == 0].sum()),
        "precision": precision,
        "recall": recall,
        "f1": 2 * precision * recall / (precision + recall) if precision + recall else 0.0,
        "null_weight": float(null.sum()),
        "known_weight": float(known.sum()),
        "null_absolute_error": float(abs(null.sum() - null_truth)),
        "cosine": float(np.dot(reconstructed, observed) / norm) if norm else 0.0,
        "residual_l1": float(np.abs(reconstructed - observed).sum()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=25)
    parser.add_argument("--replicates", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20261023)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--cancers", nargs="+", default=list(CANCERS))
    parser.add_argument("--spa", action="store_true")
    args = parser.parse_args()
    args.source = args.source.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    ref, gen, counts, truth, metadata = generate(
        args.source, args.samples, args.replicates, args.seed, args.cancers
    )
    matrix, names = ref.to_numpy(), ref.columns.tolist()
    metadata.to_csv(args.output / "samples.tsv", sep="\t", index=False)
    pd.DataFrame(counts.T, index=ref.index, columns=metadata["sample"]).to_csv(
        args.output / "catalogs.tsv", sep="\t"
    )
    ref.to_csv(args.output / "reference.tsv", sep="\t")
    np.savez_compressed(
        args.output / "inputs.npz",
        counts=counts,
        truth=truth,
        reference=matrix,
        generating=gen.to_numpy(),
        signature_names=np.array(names),
        generating_names=gen.columns.to_numpy(dtype=str),
    )
    provenance = {
        "source_url": "https://github.com/8medom/SigFitTest",
        "source_commit": subprocess.check_output(
            ["git", "-C", str(args.source), "rev-parse", "HEAD"], text=True
        ).strip(),
        "mutagene_base_commit": subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True
        ).strip(),
        "fitter_sha256": sha256(Path(__file__).with_name("frozen_identify.py")),
        "runner_sha256": sha256(__file__),
        "seed": args.seed,
        "samples_per_cell": args.samples,
        "replicates": args.replicates,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "versions": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "pandas")},
        "input_hashes": {
            str(p.relative_to(args.source)): sha256(p)
            for p in [
                args.source / "input/COSMIC_v3_SBS_GRCh38.txt",
                args.source / "input/COSMIC_v3.3.1_SBS_GRCh38.txt",
                *[
                    args.source / f"cosmic tissue data/signature_contributions_WGS_{c}.dat"
                    for c in args.cancers
                ],
            ]
        },
        "notes": [
            "Reduced replication of the 2024 design, not an exact full-paper reproduction.",
            "Uniform and substitution-specific noise are additional controls.",
            "Scaled MLE is an experimental objective-rescaling ablation.",
            "No exposure renormalization follows fitting or the ten-mutation reporting cutoff.",
        ],
    }
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    rows, exposure_rows, null_rows = [], [], []
    tasks = ((i, c, matrix, names) for i, c in enumerate(counts))
    start = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")
    ) as pool:
        for completed, result in enumerate(pool.map(fit_one, tasks, chunksize=4), start=1):
            for index, method, known, null, elapsed, diagnostics in result:
                row = metadata.iloc[index].to_dict()
                row.update(method=method, seconds=elapsed, **diagnostics)
                row.update(
                    score(
                        known,
                        null,
                        truth[index],
                        counts[index],
                        names,
                        gen.columns.tolist(),
                        matrix,
                        row["null_truth"],
                    )
                )
                rows.append(row)
                exposure_rows.append(known)
                null_rows.append(null)
            if completed % 50 == 0 or completed == len(counts):
                print(
                    f"MutaGene: {completed}/{len(counts)} samples, {time.perf_counter() - start:.1f}s",
                    flush=True,
                )
                pd.DataFrame(rows).to_csv(args.output / "metrics.tsv", sep="\t", index=False)
    if args.spa:
        from SigProfilerAssignment import Analyzer
        from sigProfilerPlotting import get_context_reference

        # SPA reorders samples, but its custom reference keeps the supplied row order.
        spa_reference = args.output / "spa_reference.tsv"
        ref.loc[get_context_reference("96")].to_csv(
            spa_reference, sep="\t", index_label="MutationType"
        )
        spa_start = time.perf_counter()
        Analyzer.cosmic_fit(
            samples=str(args.output / "catalogs.tsv"),
            output=str(args.output / "spa"),
            signature_database=str(spa_reference),
            input_type="matrix",
            genome_build="GRCh38",
            cosmic_version=3,
            make_plots=False,
            sample_reconstruction_plots=False,
            export_probabilities=False,
            cpu=args.workers,
        )
        spa_seconds = time.perf_counter() - spa_start
        activities = pd.read_csv(
            args.output / "spa/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
            sep="\t",
            index_col=0,
        )
        for index, row in metadata.iterrows():
            known = (
                activities.loc[row["sample"]].reindex(names, fill_value=0).to_numpy(dtype=float)
                / row["burden"]
            )
            null = np.zeros(6)
            result = row.to_dict()
            result.update(method="sigprofilerassignment", seconds=spa_seconds / len(metadata))
            result.update(
                score(
                    known,
                    null,
                    truth[index],
                    counts[index],
                    names,
                    gen.columns.tolist(),
                    matrix,
                    row["null_truth"],
                )
            )
            rows.append(result)
            exposure_rows.append(known)
            null_rows.append(null)
        provenance["versions"]["SigProfilerAssignment"] = importlib.metadata.version(
            "SigProfilerAssignment"
        )
    metrics = pd.DataFrame(rows)
    metrics.to_csv(args.output / "metrics.tsv", sep="\t", index=False)
    np.savez_compressed(
        args.output / "estimates.npz", known=np.array(exposure_rows), null=np.array(null_rows)
    )
    metrics.groupby(["condition", "burden", "method"]).mean(numeric_only=True).to_csv(
        args.output / "summary.tsv", sep="\t"
    )
    provenance["elapsed_seconds"] = time.perf_counter() - start
    provenance["output_hashes"] = {
        p.name: sha256(p)
        for p in args.output.iterdir()
        if p.is_file() and p.name != "provenance.json"
    }
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Finished {len(metrics)} fits", flush=True)


if __name__ == "__main__":
    main()
