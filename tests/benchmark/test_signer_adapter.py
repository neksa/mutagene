"""Protect signature alignment and native-to-relative exposure conversion."""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SOURCE = Path(__file__).resolve().parents[2] / "analysis/signature-benchmark-2026/signer_adapter.py"
spec = importlib.util.spec_from_file_location("signer_adapter", SOURCE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixture(tmp_path, labels=("A", "B", "background")):
    prefix = tmp_path / "0001"
    info = dict(
        error=None,
        sample="sample",
        task_id=1,
        background=True,
        weights=[0.2, 0.3, 0.5],
        raw_total=50,
    )
    prefix.with_suffix(".json").write_text(json.dumps(info))
    pd.DataFrame(dict(signature=labels, median_count=[10, 15, 25], weight=[0.2, 0.3, 0.5])).to_csv(
        prefix.with_suffix(".estimates.tsv"), sep="\t", index=False
    )
    task = dict(task_id=1, sample="sample", background=True)
    return prefix, task


def test_relative_weights_preserve_native_exposure_ratios(tmp_path):
    prefix, task = fixture(tmp_path)
    weights, native, info = module.read_fit(prefix, task, ["A", "B"], ["background"])
    np.testing.assert_allclose(weights, [0.2, 0.3, 0.5])
    np.testing.assert_allclose(native, [10, 15, 25])
    assert info["raw_total"] == 50


@pytest.mark.parametrize("labels", [("B", "A", "background"), ("A", "A", "background")])
def test_misaligned_signature_labels_fail(tmp_path, labels):
    prefix, task = fixture(tmp_path, labels)
    with pytest.raises(ValueError, match="labels or column order"):
        module.read_fit(prefix, task, ["A", "B"], ["background"])


def test_task_mismatch_fails(tmp_path):
    prefix, task = fixture(tmp_path)
    task["sample"] = "another"
    with pytest.raises(ValueError, match="identifiers"):
        module.read_fit(prefix, task, ["A", "B"], ["background"])


def test_sampler_error_does_not_create_an_estimate(tmp_path):
    prefix, task = fixture(tmp_path)
    prefix.with_suffix(".json").write_text(json.dumps(dict(error="Sampler failure")))
    with pytest.raises(RuntimeError, match="Sampler failure"):
        module.read_fit(prefix, task, ["A", "B"], ["background"])
