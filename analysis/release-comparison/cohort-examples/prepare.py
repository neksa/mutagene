"""Prepare labelled PCAWG counts or TOPACIO event counts from the public source files."""

import argparse
import hashlib
from collections import Counter
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
HASHES = {
    "pcawg": "90b81600a587e20fd58abd0d0396be5d08b038d620a9cf5e643573ab29528f5c",
    "topacio": "b0bbc57209e815673f1f8a8655ed22fa880347a693cfd6f9a9157da240641759",
}


def prepare(dataset, source, output):
    if hashlib.sha256(source.read_bytes()).hexdigest() != HASHES[dataset]:
        raise ValueError("The source checksum differs from the study input")
    if dataset == "pcawg":
        counts = pd.read_csv(source, comment="#", sep="\t", index_col=0)
        channels = pd.read_csv(
            HERE.parent / "performance/reference.tsv", sep="\t", index_col=0
        ).index
        if set(counts.index) != set(channels) or counts.shape != (96, 146):
            raise ValueError("Expected 146 labelled PCAWG profiles")
        counts.loc[channels].to_csv(output, sep="\t")
        return

    from openpyxl import load_workbook

    sheet = load_workbook(source, read_only=True, data_only=True)["Figure 1"]
    rows = iter(sheet.values)
    columns = next(rows)
    counts = Counter()
    for row in rows:
        patient = dict(zip(columns, row))
        status = patient["SigMA"]
        if status is None:
            continue
        group = "Positive" if status else "Negative"
        counts[group, patient["PFS"], patient["PFS_Censor"]] += 1
    events = pd.DataFrame(
        [(group, day, count, censor) for (group, day, censor), count in counts.items()],
        columns=["SigMA", "days", "count", "censor"],
    )
    events.sort_values(["SigMA", "days", "censor"]).to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", choices=HASHES)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    prepare(args.dataset, args.source, args.output)
