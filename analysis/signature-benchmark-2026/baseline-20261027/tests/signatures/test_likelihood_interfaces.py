"""Verify shared MLE behavior through the CLI and web application."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from mutagene.cli.signature_menu import SignatureMenu
from mutagene.io.matrix import write_count_matrix
from mutagene.io.profile import read_signatures
from mutagene.signatures.identify import decompose_mutational_profile_counts
from mutagene.webapp.analysis import profile_dict_to_array, run_cohort_analysis
from tests.cli.local_fixtures import LOCAL_GENOME, LOCAL_MAF


def test_cli_default_matches_direct_python_likelihood(tmp_path):
    counts = np.arange(1, 97)
    path = tmp_path / "matrix.tsv"
    with path.open("w") as stream:
        write_count_matrix(stream, {"sample": counts})
    output = tmp_path / "diagnostics.json"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "mutagene",
            "signature",
            "-f",
            "MATRIX",
            "-i",
            str(path),
            "-s",
            "COSMICv3",
            "-o",
            str(tmp_path / "fit.tsv"),
            "--diagnostics-out",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    data = json.loads(output.read_text())
    direct = {}
    weights, _, _ = decompose_mutational_profile_counts(
        counts, read_signatures("COSMICv3"), diagnostics=direct
    )
    assert data["method"] == "MLE"
    assert data["samples"]["sample"]["nll_per_mutation"] == pytest.approx(
        direct["nll_per_mutation"], abs=1e-10
    )
    assert data["samples"]["sample"]["null_components"] == pytest.approx(weights[-6:], abs=1e-8)


@pytest.fixture
def local_web_genome(monkeypatch):
    monkeypatch.setattr(
        "mutagene.webapp.analysis.GenomeManager.get_genome_path",
        lambda self, genome: Path(LOCAL_GENOME),
    )


def test_web_pooled_and_per_sample_fits_use_certified_mle(tmp_path, local_web_genome):
    result = run_cohort_analysis(Path(LOCAL_MAF), tmp_path / "web")
    assert result["analysis_status"]["signatures"] == "complete"
    summary = result["signatures"]
    assert summary["method"] == "MLE"
    info = {}
    decompose_mutational_profile_counts(
        profile_dict_to_array(result["profiles"]), read_signatures("COSMICv3"), diagnostics=info
    )
    assert summary["diagnostics"]["nll_per_mutation"] == pytest.approx(
        info["nll_per_mutation"], abs=1e-10
    )
    assert len(summary["per_sample"]) == 3
    for sample in summary["per_sample"].values():
        assert sample["diagnostics"]["estimator"] == "multinomial_mle"
        assert sample["diagnostics"]["fallback"] is None
        assert sample["diagnostics"]["optimizer_success"]


def test_web_marks_uncertified_fit_as_failed(tmp_path, local_web_genome, monkeypatch):
    def failed(*args, **kwargs):
        raise RuntimeError("Likelihood optimization fails certification")

    monkeypatch.setattr("mutagene.signatures.identify.fit_likelihood", failed)
    result = run_cohort_analysis(Path(LOCAL_MAF), tmp_path / "web")
    assert result["analysis_status"]["signatures"] == "failed"
    assert "certification" in result["signatures"]["error"]


def test_cli_returns_failure_for_uncertified_fit(tmp_path, monkeypatch):
    path = tmp_path / "matrix.tsv"
    with path.open("w") as stream:
        write_count_matrix(stream, {"sample": np.ones(96)})
    parser = argparse.ArgumentParser()
    menu = SignatureMenu(parser)
    args = parser.parse_args(["-f", "MATRIX", "-i", str(path), "-o", str(tmp_path / "fit.tsv")])

    def failed(*args, **kwargs):
        raise RuntimeError("Likelihood optimization fails certification")

    monkeypatch.setattr("mutagene.signatures.identify.fit_likelihood", failed)
    try:
        with pytest.raises(SystemExit) as error:
            menu.identify(args)
        assert error.value.code == 1
    finally:
        args.infile.close()
        args.outfile.close()
