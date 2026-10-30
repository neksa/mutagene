"""Validate native signeR outputs and retain explicit relative-exposure semantics."""

import json
from pathlib import Path

import numpy as np
import pandas as pd


def read_fit(prefix, task, known_names, background_names):
    """Return normalized weights, native median counts, and sampler diagnostics."""
    prefix = Path(prefix)
    info = json.loads(prefix.with_suffix(".json").read_text())
    if info["error"] is not None:
        raise RuntimeError(info["error"])
    if info["sample"] != task["sample"] or info["task_id"] != int(task["task_id"]):
        raise ValueError("signeR task identifiers do not match")
    background = str(info["background"]).lower() == "true"
    if background != bool(task["background"]):
        raise ValueError("signeR background configuration does not match")
    expected = list(known_names) + (list(background_names) if background else [])
    estimates = pd.read_csv(prefix.with_suffix(".estimates.tsv"), sep="\t")
    if estimates.signature.tolist() != expected:
        raise ValueError("signeR signature labels or column order differ")
    native = estimates.median_count.to_numpy(float)
    weights = estimates.weight.to_numpy(float)
    if not np.isfinite(native).all() or (native < 0).any() or native.sum() <= 0:
        raise ValueError("signeR median exposures must be finite and nonnegative")
    if not np.allclose(weights, native / native.sum(), atol=1e-12, rtol=0):
        raise ValueError("signeR relative exposures do not match native median normalization")
    if not np.allclose(weights, info["weights"], atol=1e-12, rtol=0):
        raise ValueError("signeR estimate records disagree")
    if not np.isclose(native.sum(), info["raw_total"], atol=1e-9, rtol=1e-12):
        raise ValueError("signeR native exposure totals disagree")
    return weights, native, info
