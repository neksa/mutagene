"""Inspect optimizer candidates that the production MLE fallback discards."""

import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from audit_likelihood import OUT, ROOT, certificate
from run_benchmark import sha256

from mutagene.signatures import get_dummy_signatures_lists, identify


def inspect(task):
    sample, burden, background, counts, reference, names = task
    original = identify.minimize
    captured = []

    def capture(*args, **kwargs):
        result = original(*args, **kwargs)
        captured.append(result)
        return result

    identify.minimize = capture
    info = {}
    try:
        returned, _, _ = identify.decompose_mutational_profile_counts(
            counts, (reference, names), func="MLE", enable_dummy=background, diagnostics=info
        )
    finally:
        identify.minimize = original
    candidate = captured[0].x
    dummy = np.array([values for _, values in get_dummy_signatures_lists()]).T
    matrix = np.column_stack([reference, dummy]) if background else reference
    candidate_nll, candidate_gap = certificate(matrix, counts, candidate)
    returned_nll, returned_gap = certificate(matrix, counts, returned)
    return {
        "sample": sample,
        "burden": burden,
        "background": background,
        "optimizer_success": info["optimizer_success"],
        "candidate_weight_sum": candidate.sum(),
        "candidate_min_weight": candidate.min(),
        "candidate_nll_per_mutation": candidate_nll,
        "candidate_gap_bound": candidate_gap,
        "returned_nll_per_mutation": returned_nll,
        "returned_gap_bound": returned_gap,
        "fallback_nll_increase": returned_nll - candidate_nll,
    }


def main():
    source = ROOT / "sparse-results/evaluate"
    spec = json.loads((source / "specification.json").read_text())
    assert (
        sha256(ROOT.parents[1] / "mutagene/signatures/identify.py")
        == spec["source_sha256"]["identify.py"]
    )
    metrics = pd.read_csv(source / "metrics.tsv", sep="\t", low_memory=False)
    samples = pd.read_csv(source / "samples.tsv", sep="\t")
    indices = dict(zip(samples["sample"], samples.index))
    with np.load(source / "inputs.npz") as archive:
        data = {key: archive[key] for key in archive.files}
    failures = metrics[(metrics.method == "MLE") & metrics.fallback.notna()]
    tasks = [
        (
            r.sample,
            r.burden,
            r.background,
            data["counts"][indices[r.sample]],
            data["reference"],
            data["names"].tolist(),
        )
        for r in failures.itertuples()
    ]
    with ProcessPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(inspect, tasks))
    pd.DataFrame(rows).to_csv(OUT / "discarded-candidates.tsv", sep="\t", index=False)
    print(f"Inspected {len(rows)} archived MLE failures")


if __name__ == "__main__":
    main()
