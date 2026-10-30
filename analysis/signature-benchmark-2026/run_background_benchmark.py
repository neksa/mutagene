"""Reassess repaired fitting on paired profiles and a prespecified fresh simulation panel."""

import argparse
import importlib.metadata
import importlib.util
import json
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from run_benchmark import CANCERS, generate, score, sha256

from mutagene.signatures.identify import NegLogLik, decompose_mutational_profile_counts

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "background-results"


def fit_sample(task):
    index, counts, matrix, names, configs = task
    rows = []
    for label, method, background, policy in configs:
        info = {}
        begin = time.perf_counter()
        try:
            weights, _, _ = decompose_mutational_profile_counts(
                counts,
                (matrix, names),
                func=method,
                enable_dummy=background,
                others_threshold=0,
                diagnostics=info,
                background_selection=policy,
            )
            known = weights[: len(names)]
            null = weights[len(names) :] if background else np.zeros(6)
            error = None
        except Exception as exc:
            known, null, error = None, None, repr(exc)
        rows.append(
            (index, label, background, known, null, time.perf_counter() - begin, info, error)
        )
    return rows


CONFIGS = [
    ("joint-BIC", "sparse-bic", True, "joint"),
    ("block-BIC", "sparse-bic", True, "block"),
    ("fixed-BIC", "sparse-bic", True, "fixed"),
    ("no-background-BIC", "sparse-bic", False, "fixed"),
    ("MLE", "MLE", True, "fixed"),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["paired", "fresh"])
    parser.add_argument("--source", type=Path, default=Path("/tmp/mutagene-sigfittest-2024"))
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    out = OUT / args.phase
    if out.exists():
        raise ValueError("Choose a fresh output directory; never overwrite benchmark evidence")
    commit = subprocess.check_output(
        ["git", "-C", str(args.source), "rev-parse", "HEAD"], text=True
    ).strip()
    assert commit == "e94bac6e9b04fb13500d1904b5e0e8d58b00d8ad"
    upstream_hashes = json.loads((ROOT / "results/provenance.json").read_text())["input_hashes"]
    for filename, digest in upstream_hashes.items():
        assert sha256(args.source / filename) == digest
    old = ROOT / "repaired-results/fresh"
    if args.phase == "paired":
        with np.load(old / "inputs.npz") as archive:
            inputs = {name: archive[name] for name in archive.files}
        metadata = pd.read_csv(old / "samples.tsv", sep="\t")
        from run_benchmark import channel_order

        reference = pd.DataFrame(
            inputs["reference"], index=channel_order(), columns=inputs["names"]
        )
        counts, truth, generating_names = (
            inputs["counts"],
            inputs["truth"],
            inputs["generating_names"].tolist(),
        )
        seed = 20261027
    else:
        seed = 202610273
        reference, generating, counts, truth, metadata = generate(args.source, 8, 1, seed, CANCERS)
        metadata["sample"] = "fresh_" + metadata["sample"]
        generating_names = generating.columns.tolist()
    out.mkdir(parents=True)
    metadata.to_csv(out / "samples.tsv", sep="\t", index=False)
    reference.to_csv(out / "reference.tsv", sep="\t")
    pd.DataFrame(counts.T, index=reference.index, columns=metadata["sample"]).to_csv(
        out / "catalogs.tsv", sep="\t", index_label="MutationType"
    )
    np.savez_compressed(
        out / "inputs.npz",
        counts=counts,
        truth=truth,
        reference=reference.to_numpy(),
        names=reference.columns.to_numpy(dtype=str),
        generating_names=np.array(generating_names, dtype=str),
    )
    specification = {
        "phase": args.phase,
        "seed": seed,
        "profiles": len(counts),
        "methods": CONFIGS + ["SigProfilerAssignment"],
        "source_commit": commit,
        "upstream_input_hashes": upstream_hashes,
        "source_sha256": {
            name: sha256(ROOT.parents[1] / "mutagene/signatures" / name)
            for name in ("likelihood.py", "identify.py", "sparse.py")
        },
        "baseline_source_sha256": sha256(
            ROOT / "baseline-20261026/mutagene/signatures/identify.py"
        ),
        "runner_sha256": sha256(__file__),
        "versions": {
            name: importlib.metadata.version(name)
            for name in ("numpy", "scipy", "pandas", "SigProfilerAssignment")
        },
        "frozen_sparse_criterion": json.loads((ROOT / "sparse-results/selection.json").read_text()),
        "prespecified_assessment": "No retuning. Compare joint and block background BIC against fixed BIC, no-background BIC, MLE, and SPA. Report attribution error, background mass, all conditions and burdens. Select a default relative to MLE; disclose low-burden and matched-noise limitations.",
        "limitations": "Paired profiles inform the repair. Fresh draws share empirical exposure distributions and are not independent patient data.",
    }
    (out / "specification.json").write_text(json.dumps(specification, indent=2) + "\n")
    rows, knowns, nulls = [], [], []

    def record(group):
        for index, method, background, known, null, seconds, info, error in group:
            row = metadata.iloc[index].to_dict()
            row.update(method=method, background=background, seconds=seconds, error=error)
            row.update(info)
            if error is None:
                row.update(
                    score(
                        known,
                        null,
                        truth[index],
                        counts[index],
                        reference.columns.tolist(),
                        generating_names,
                        reference.to_numpy(),
                        row["null_truth"],
                    )
                )
                row["active_known"] = int(np.count_nonzero(known * counts[index].sum() >= 10))
                row["positive_known"] = int(np.count_nonzero(known > 1e-7))
                from mutagene.signatures import get_dummy_signatures_lists

                dummy = np.array([values for _, values in get_dummy_signatures_lists()]).T
                full = np.column_stack([reference.to_numpy(), dummy])
                row["count_log_likelihood"] = -NegLogLik(np.r_[known, null], full, counts[index])
            else:
                known, null = np.full(reference.shape[1], np.nan), np.full(6, np.nan)
            rows.append(row)
            knowns.append(known)
            nulls.append(null)

    begin = time.perf_counter()
    tasks = [
        (i, counts[i], reference.to_numpy(), reference.columns.tolist(), CONFIGS)
        for i in range(len(counts))
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for completed, group in enumerate(pool.map(fit_sample, tasks), 1):
            record(group)
            if completed % 32 == 0:
                pd.DataFrame(rows).to_csv(out / "metrics.tsv", sep="\t", index=False)
                print(
                    f"{args.phase}: {completed}/{len(counts)} profiles, {time.perf_counter()-begin:.1f}s",
                    flush=True,
                )
    from SigProfilerAssignment import Analyzer
    from sigProfilerPlotting import get_context_reference

    reference.loc[get_context_reference("96")].to_csv(
        out / "spa_reference.tsv", sep="\t", index_label="MutationType"
    )
    start = time.perf_counter()
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
    elapsed = time.perf_counter() - start
    activities = pd.read_csv(
        out / "spa/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
        sep="\t",
        index_col=0,
    )
    for i in range(len(counts)):
        known = (
            activities.loc[metadata.iloc[i]["sample"]].reindex(reference.columns).to_numpy(float)
            / counts[i].sum()
        )
        record(
            [
                (
                    i,
                    "SigProfilerAssignment",
                    False,
                    known,
                    np.zeros(6),
                    elapsed / len(counts),
                    {},
                    None,
                )
            ]
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "metrics.tsv", sep="\t", index=False)
    np.savez_compressed(out / "estimates.npz", known=np.array(knowns), null=np.array(nulls))
    frame.groupby(["condition", "burden", "method", "background"]).mean(numeric_only=True).to_csv(
        out / "summary.tsv", sep="\t"
    )
    print(
        f"Finished {args.phase}: {len(frame)} fits, {frame.error.notna().sum()} errors", flush=True
    )


if __name__ == "__main__":
    main()
