"""Check installed source behavior against every frozen validation estimate."""

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1]))

from mutagene.signatures import get_dummy_signatures_lists  # noqa: E402
from mutagene.signatures.sparse import fit_sparse  # noqa: E402


def verify(task):
    matrix, counts, known_count, expected = task
    info = {}
    actual = fit_sparse(matrix, counts, known_count, diagnostics=info)
    np.testing.assert_array_equal(actual, expected)
    assert info["maximum_refit_optimality_gap"] <= 2e-5
    assert info["fallback"] is None
    return info["search_start"], info["maximum_refit_optimality_gap"]


def main():
    phase_records = {}
    for phase in ["validation", "sensitivity"]:
        folder = ROOT / phase
        data = dict(np.load(folder / "inputs.npz"))
        meta = pd.read_csv(folder / "samples.tsv", sep="\t")
        saved = np.load(folder / "estimates.npz")["weights"].reshape(len(meta), 2, -1)
        dummy = np.array([v for _, v in get_dummy_signatures_lists()]).T
        matrix = np.column_stack([data["reference"], dummy])
        tasks = [
            (matrix, counts, len(data["names"]), saved[i, 1])
            for i, counts in enumerate(data["counts"])
        ]
        records = []
        with ProcessPoolExecutor(max_workers=4) as pool:
            for done, result in enumerate(pool.map(verify, tasks), 1):
                records.append(result)
                if done % 64 == 0:
                    print(f"{phase}: verified {done}/{len(tasks)} production estimates", flush=True)
        phase_records[phase] = dict(
            profiles=len(records),
            exact_estimate_matches=len(records),
            full_starts=sum(r[0] == "full" for r in records),
            single_starts=sum(r[0] == "single" for r in records),
            maximum_refit_optimality_gap=max(r[1] for r in records),
        )
    (ROOT / "results/implementation-verification.json").write_text(
        json.dumps(phase_records, indent=2) + "\n"
    )
    print(json.dumps(phase_records, indent=2), flush=True)


if __name__ == "__main__":
    main()
