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
