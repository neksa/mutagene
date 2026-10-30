"""Tune sparse selection, freeze its criterion, and evaluate on fresh synthetic profiles."""

import argparse
import importlib.metadata
import json
import multiprocessing
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from run_benchmark import CANCERS, generate, score, sha256

from mutagene.signatures import identify

ROOT = Path(__file__).resolve().parent
LEGACY = (
    "MLE",
    "MLE-scaled",
    "Frobenius",
    "FrobeniusZero",
    "Cos",
    "Elastic",
    "KL",
    "JS",
    "AICc",
    "BIC",
    "Compat",
)
SPARSE = ("sparse-aicc", "sparse-bic")


def fit_sample(task):
    index, counts, matrix, names, methods = task
    rows = []
    for method in methods:
        for background in (True, False):
            info = {}
            start = time.perf_counter()
            try:
                weights, _, _ = identify.decompose_mutational_profile_counts(
                    counts,
                    (matrix, names),
                    func=method,
                    enable_dummy=background,
                    others_threshold=0,
                    diagnostics=info,
                )
                known = weights[: len(names)]
                null = weights[len(names) :] if background else np.zeros(6)
                rows.append(
                    (
                        index,
                        method,
                        background,
                        known,
                        null,
                        time.perf_counter() - start,
                        info,
                        None,
                    )
                )
            except Exception as error:
                rows.append(
                    (
                        index,
                        method,
                        background,
                        None,
                        None,
                        time.perf_counter() - start,
                        info,
                        str(error),
                    )
                )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["tune", "evaluate"])
    parser.add_argument("--source", type=Path, default=Path("/tmp/mutagene-sigfittest-2024"))
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    phase = args.phase
    out = ROOT / "sparse-results" / phase
    if out.exists():
        raise ValueError("Choose a fresh phase directory; this runner does not overwrite evidence")
    if phase == "evaluate" and not (ROOT / "sparse-results/selection.json").exists():
        raise ValueError("Freeze the tuning selection before evaluation")
    seed = 202610231 if phase == "tune" else 202610232
    samples = 2 if phase == "tune" else 8
    methods = SPARSE if phase == "tune" else (*LEGACY, *SPARSE)
    ref, gen, counts, truth, metadata = generate(args.source, samples, 1, seed, CANCERS)
    metadata["sample"] = phase + "_" + metadata["sample"]
    out.mkdir(parents=True)
    metadata.to_csv(out / "samples.tsv", sep="\t", index=False)
    ref.to_csv(out / "reference.tsv", sep="\t")
    pd.DataFrame(counts.T, index=ref.index, columns=metadata["sample"]).to_csv(
        out / "catalogs.tsv", sep="\t", index_label="MutationType"
    )
    np.savez_compressed(
        out / "inputs.npz",
        counts=counts,
        truth=truth,
        reference=ref.to_numpy(),
        names=ref.columns.to_numpy(dtype=str),
        generating_names=gen.columns.to_numpy(dtype=str),
    )
    specification = {
        "phase": phase,
        "seed": seed,
        "samples_per_cancer_burden_condition": samples,
        "profiles": len(counts),
        "methods": list(methods),
        "background": [True, False],
        "workers": args.workers,
        "source_commit": "e94bac6e9b04fb13500d1904b5e0e8d58b00d8ad",
        "source_sha256": {
            name: sha256(ROOT.parents[1] / "mutagene/signatures" / name)
            for name in ["sparse.py", "identify.py"]
        },
        "runner_sha256": sha256(__file__),
        "versions": {
            p: importlib.metadata.version(p)
            for p in ["numpy", "scipy", "pandas", "SigProfilerAssignment"]
        },
        "selection_rule": "Choose the lowest mean aggregate-null TAE across both background policies and all tuning cells; require zero failed fits.",
        "limitations": "Fresh Monte Carlo draws share empirical exposure distributions; this is not independent patient validation.",
    }
    if phase == "evaluate":
        specification["frozen_selection"] = json.loads(
            (ROOT / "sparse-results/selection.json").read_text()
        )
    (out / "specification.json").write_text(json.dumps(specification, indent=2) + "\n")
    rows, knowns, nulls = [], [], []
    tasks = ((i, c, ref.to_numpy(), ref.columns.tolist(), methods) for i, c in enumerate(counts))
    start = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")
    ) as pool:
        for completed, group in enumerate(pool.map(fit_sample, tasks, chunksize=1), 1):
            for index, method, background, known, null, elapsed, info, error in group:
                row = metadata.iloc[index].to_dict()
                trace = info.pop("selection_trace", [])
                info.pop("selected_known_indices", None)
                row.update(
                    method=method, background=background, seconds=elapsed, error=error, **info
                )
                row["selection_steps"] = len(trace)
                if error is None:
                    row.update(
                        score(
                            known,
                            null,
                            truth[index],
                            counts[index],
                            ref.columns.tolist(),
                            gen.columns.tolist(),
                            ref.to_numpy(),
                            row["null_truth"],
                        )
                    )
                    row["active_known"] = int(np.count_nonzero(known * counts[index].sum() >= 10))
                    row["positive_known"] = int(np.count_nonzero(known > 1e-7))
                else:
                    known, null = np.full(ref.shape[1], np.nan), np.full(6, np.nan)
                rows.append(row)
                knowns.append(known)
                nulls.append(null)
            if completed % 16 == 0 or completed == len(counts):
                pd.DataFrame(rows).to_csv(out / "metrics.tsv", sep="\t", index=False)
                print(
                    f"{phase}: {completed}/{len(counts)}, {time.perf_counter()-start:.1f}s",
                    flush=True,
                )
    if phase == "evaluate":
        from SigProfilerAssignment import Analyzer
        from sigProfilerPlotting import get_context_reference

        ref.loc[get_context_reference("96")].to_csv(
            out / "spa_reference.tsv", sep="\t", index_label="MutationType"
        )
        begin = time.perf_counter()
        Analyzer.cosmic_fit(
            samples=str(out / "catalogs.tsv"),
            output=str(out / "spa"),
            signature_database=str(out / "spa_reference.tsv"),
            input_type="matrix",
            genome_build="GRCh38",
            cosmic_version=3,
            make_plots=False,
            sample_reconstruction_plots=False,
            export_probabilities=False,
            cpu=args.workers,
        )
        elapsed = time.perf_counter() - begin
        activities = pd.read_csv(
            out / "spa/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
            sep="\t",
            index_col=0,
        )
        for i, row in metadata.iterrows():
            known = (
                activities.loc[row["sample"]].reindex(ref.columns).to_numpy(float) / row["burden"]
            )
            null = np.zeros(6)
            result = row.to_dict()
            result.update(
                method="SigProfilerAssignment",
                background=False,
                seconds=elapsed / len(metadata),
                error=None,
            )
            result.update(
                score(
                    known,
                    null,
                    truth[i],
                    counts[i],
                    ref.columns.tolist(),
                    gen.columns.tolist(),
                    ref.to_numpy(),
                    row["null_truth"],
                )
            )
            result.update(
                active_known=int(np.count_nonzero(known * row["burden"] >= 10)),
                positive_known=int(np.count_nonzero(known > 1e-7)),
            )
            rows.append(result)
            knowns.append(known)
            nulls.append(null)
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "metrics.tsv", sep="\t", index=False)
    np.savez_compressed(out / "estimates.npz", known=np.array(knowns), null=np.array(nulls))
    frame.groupby(["condition", "burden", "method", "background"]).mean(numeric_only=True).to_csv(
        out / "summary.tsv", sep="\t"
    )
    if phase == "tune":
        candidates = {}
        for method, group in frame.groupby("method"):
            candidates[method] = {
                "mean_error": float(group.aggregate_null_tae.mean()),
                "failures": int(group.error.notna().sum()),
            }
        eligible = [m for m in candidates if not candidates[m]["failures"]]
        if not eligible:
            raise RuntimeError("No sparse criterion completes tuning without failures")
        selected = min(eligible, key=lambda m: candidates[m]["mean_error"])
        (ROOT / "sparse-results/selection.json").write_text(
            json.dumps(
                {
                    "selected": selected,
                    "candidates": candidates,
                    "tuning_metrics_sha256": sha256(out / "metrics.tsv"),
                    "sparse_source_sha256": specification["source_sha256"]["sparse.py"],
                },
                indent=2,
            )
            + "\n"
        )
    print(f"Finished {phase}: {len(frame)} fits", flush=True)


if __name__ == "__main__":
    main()
