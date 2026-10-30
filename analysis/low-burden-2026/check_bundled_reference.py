"""Check the unchanged policy with the package's bundled COSMICv3 reference."""

import importlib.util
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1]))
sys.path.insert(0, str(ROOT.parent / "signature-benchmark-2026"))

from run_benchmark import score  # noqa: E402

from mutagene.io.profile import read_signatures  # noqa: E402
from mutagene.signatures import get_dummy_signatures_lists  # noqa: E402
from mutagene.signatures.sparse import fit_sparse  # noqa: E402

spec = importlib.util.spec_from_file_location(
    "previous_sparse", ROOT / "baseline-source/mutagene/signatures/sparse.py"
)
previous = importlib.util.module_from_spec(spec)
spec.loader.exec_module(previous)


def evaluate(task):
    matrix, counts, truth, names, generating, metadata = task
    rows = []
    for name, fit in [("joint-BIC", previous.fit_sparse), ("adaptive-BIC", fit_sparse)]:
        weights = fit(matrix, counts, len(names))
        metrics = score(
            weights[: len(names)],
            weights[len(names) :],
            truth,
            counts,
            names,
            generating,
            matrix[:, : len(names)],
            metadata["null_truth"],
        )
        rows.append(dict(metadata, method=name, **metrics))
    return rows


def main():
    out = ROOT / "bundled-reference"
    out.mkdir()
    data = dict(np.load(ROOT / "validation/inputs.npz"))
    meta = pd.read_csv(ROOT / "validation/samples.tsv", sep="\t")
    indices = meta.index[(meta.sample_index == 0) & (meta.burden <= 1000)]
    reference, names = read_signatures("COSMICv3")
    names = [name.removeprefix("COSMICv3-") for name in names]
    assert names == data["names"].tolist()
    matrix = np.column_stack([reference, np.array([v for _, v in get_dummy_signatures_lists()]).T])
    matrix /= matrix.sum(axis=0)
    (out / "specification.json").write_text(
        json.dumps(
            dict(
                profiles=len(indices),
                selection="sample_index == 0 across every cancer, condition, and low burden; no further tuning",
                reference="Bundled COSMICv3",
            ),
            indent=2,
        )
        + "\n"
    )
    tasks = [
        (
            matrix,
            data["counts"][i],
            data["truth"][i],
            names,
            data["generating_names"].tolist(),
            meta.iloc[i].to_dict(),
        )
        for i in indices
    ]
    rows = []
    with ProcessPoolExecutor(max_workers=2) as pool:
        for group in pool.map(evaluate, tasks):
            rows.extend(group)
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "metrics.tsv", sep="\t", index=False)
    summary = frame.groupby(["burden", "method"])[["aggregate_null_tae", "paper_tae"]].mean()
    summary.to_csv(out / "summary.tsv", sep="\t")
    print(summary.to_string())


if __name__ == "__main__":
    main()
