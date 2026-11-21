"""Render the original benchmark and cohort figures from their numerical inputs."""

import argparse
import runpy
from pathlib import Path

import pandas as pd

PAPER = Path(__file__).resolve().parents[1]
BASE = PAPER.parent / "analysis/release-comparison"
OUT = PAPER / "figures"


def render(name):
    OUT.mkdir(parents=True, exist_ok=True)
    if name == "accuracy":
        program = runpy.run_path(str(BASE / "performance/figure.py"))
        program["draw"](
            pd.read_csv(BASE / "performance/errors.tsv", sep="\t", index_col=0),
            OUT / "accuracy.pdf",
        )
    elif name == "cohort":
        program = runpy.run_path(str(BASE / "cohort-examples/figure.py"))
        program["draw"](OUT / "cohort.pdf")
    elif name == "penta":
        program = runpy.run_path(str(BASE / "pcawg-penta/figure.py"))
        counts = pd.read_csv(BASE / "pcawg-penta/pcawg-counts.tsv.gz", sep="\t", index_col=0)
        program["waterfall"](counts, OUT / "penta.pdf")
    elif name == "contexts":
        program = runpy.run_path(str(BASE / "extended-contexts/example.py"))
        p = pd.read_csv(PAPER / "inputs/SBS1536-display.tsv", sep="\t", index_col=0).SBS39_P
        indel = pd.read_csv(PAPER / "inputs/ID83-display.tsv", sep="\t", index_col=0).ID1
        program["draw"](p / p.sum(), indel / indel.sum(), OUT)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name", choices=["accuracy", "cohort", "penta", "contexts"])
    render(parser.parse_args().name)
