"""Report signature command failures without success diagnostics."""

import argparse

import numpy as np
import pytest

from mutagene.cli.signature_menu import SignatureMenu
from mutagene.io.matrix import write_count_matrix


def test_bootstrap_failure_returns_failure_without_success_diagnostics(tmp_path, monkeypatch):
    path = tmp_path / "matrix.tsv"
    with path.open("w") as stream:
        write_count_matrix(stream, {"sample": np.ones(96)})
    diagnostics = tmp_path / "diagnostics.json"
    parser = argparse.ArgumentParser()
    menu = SignatureMenu(parser)
    args = parser.parse_args(
        [
            "-f",
            "MATRIX",
            "-i",
            str(path),
            "--bootstrap",
            "-br",
            "10",
            "--diagnostics-out",
            str(diagnostics),
            "-o",
            str(tmp_path / "fit.tsv"),
        ]
    )
    calls = 0

    def fail_resample(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls > 1:
            raise RuntimeError("Likelihood optimization fails certification")
        return np.ones(30) / 30, None, []

    monkeypatch.setattr(
        "mutagene.cli.signature_menu.decompose_mutational_profile_counts", fail_resample
    )
    try:
        with pytest.raises(SystemExit) as error:
            menu.identify(args)
        assert error.value.code == 1
        assert not diagnostics.exists()
    finally:
        args.infile.close()
        args.outfile.close()


def test_method_option_requires_a_value():
    parser = argparse.ArgumentParser()
    SignatureMenu(parser)
    with pytest.raises(SystemExit) as error:
        parser.parse_args(["--method"])
    assert error.value.code == 2
