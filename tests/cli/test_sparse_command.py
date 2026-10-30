"""Verify the sparse command preserves background mass in its diagnostics."""

import json
import subprocess
import sys

import numpy as np

from mutagene.io.matrix import write_count_matrix


def test_uniform_profile_can_select_background_only(tmp_path):
    matrix = tmp_path / "matrix.tsv"
    with matrix.open("w") as stream:
        write_count_matrix(stream, {"uniform": np.full(96, 20)})
    diagnostics = tmp_path / "diagnostics.json"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "mutagene",
            "signature",
            "-f",
            "MATRIX",
            "-i",
            str(matrix),
            "--method",
            "sparse-bic",
            "--background-selection",
            "fixed",
            "-o",
            str(tmp_path / "fit.tsv"),
            "--diagnostics-out",
            str(diagnostics),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    info = json.loads(diagnostics.read_text())["samples"]["uniform"]
    assert info["selected_known_indices"] == []
    assert np.isclose(info["null_weight"], 1)
    assert np.isclose(info["known_signature_weight"], 0)
    assert info["fallback"] is None


def test_cli_joint_background_matches_python(tmp_path):
    from mutagene.io.profile import read_signatures
    from mutagene.signatures.identify import decompose_mutational_profile_counts

    counts = np.arange(1, 97)
    matrix = tmp_path / "matrix.tsv"
    with matrix.open("w") as stream:
        write_count_matrix(stream, {"sample": counts})
    diagnostics = tmp_path / "diagnostics.json"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "mutagene",
            "signature",
            "-f",
            "MATRIX",
            "-i",
            str(matrix),
            "-s",
            "COSMICv3",
            "--method",
            "sparse-bic",
            "--background-selection",
            "joint",
            "-o",
            str(tmp_path / "fit.tsv"),
            "--diagnostics-out",
            str(diagnostics),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    info = json.loads(diagnostics.read_text())["samples"]["sample"]
    direct = {}
    weights, _, _ = decompose_mutational_profile_counts(
        counts,
        read_signatures("COSMICv3"),
        func="sparse-bic",
        background_selection="joint",
        diagnostics=direct,
    )
    assert info["background_selection"] == "joint"
    assert info["selected_background_indices"] == direct["selected_background_indices"]
    assert info["criterion_value"] == direct["criterion_value"]
    np.testing.assert_allclose(info["null_components"], weights[-6:], atol=1e-12)
