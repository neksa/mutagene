"""Resolve immutable historical fitter sources after production repairs."""

from pathlib import Path

from run_benchmark import sha256

ROOT = Path(__file__).resolve().parent


def signature_source(name, expected):
    relative = Path("mutagene/signatures") / name
    for path in (
        ROOT.parents[1] / relative,
        ROOT / "baseline-20261027" / relative,
        ROOT / "baseline-20261026" / relative,
    ):
        if path.is_file() and sha256(path) == expected:
            return path
    raise ValueError(f"No source snapshot matches {name}: {expected}")
