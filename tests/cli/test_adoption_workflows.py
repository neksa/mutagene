"""Exercise installation-oriented workflows with offline fixtures."""

import gzip
import json
import subprocess
import sys

import pandas as pd
import pytest

from tests.cli.local_fixtures import LOCAL_GENOME, LOCAL_MAF


def run(*args):
    return subprocess.run([sys.executable, "-m", "mutagene", *args], capture_output=True, text=True)


def test_demo_and_matrix_fitting_need_no_external_genome(tmp_path):
    demo = tmp_path / "demo"
    result = run("demo", "--output-dir", str(demo))
    assert result.returncode == 0, result.stderr
    summary = json.loads((demo / "summary.json").read_text())
    assert (summary["samples"], summary["mutations"], summary["channels"]) == (3, 24, 96)
    result = run(
        "signature",
        "-f",
        "MATRIX",
        "-i",
        str(demo / "matrix.tsv"),
        "-o",
        str(tmp_path / "fit.tsv"),
        "--diagnostics-out",
        str(tmp_path / "fit.json"),
    )
    assert result.returncode == 0, result.stderr
    assert len(json.loads((tmp_path / "fit.json").read_text())["samples"]) == 3


def test_gzip_per_sample_matrix_preserves_counts(tmp_path):
    compressed = tmp_path / "input.maf.gz"
    with gzip.open(compressed, "wt") as stream, open(LOCAL_MAF) as original:
        stream.write(original.read())
    output = tmp_path / "matrix.tsv"
    result = run(
        "profile", "--per-sample", "-i", str(compressed), "-g", LOCAL_GENOME, "-o", str(output)
    )
    assert result.returncode == 0, result.stderr
    frame = pd.read_csv(output, sep="\t", index_col=0)
    assert frame.shape == (96, 3)
    assert frame.to_numpy().sum() == 24


@pytest.mark.parametrize(
    "args",
    [
        ("profile",),
        ("signature",),
        ("signature", "-i", LOCAL_MAF, "-g", LOCAL_GENOME, "--method", "invalid"),
        ("rank", "--cohort", "missing", "--cohorts-file", "/nonexistent/mutagene.tar.gz"),
    ],
)
def test_invalid_requests_return_failure(args):
    assert run(*args).returncode != 0


def test_standard_input_supports_repeated_reads(tmp_path):
    output = tmp_path / "matrix.tsv"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "mutagene",
            "profile",
            "--per-sample",
            "-i",
            "-",
            "-g",
            LOCAL_GENOME,
            "-o",
            str(output),
        ],
        input=open(LOCAL_MAF).read(),
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert pd.read_csv(output, sep="\t", index_col=0).to_numpy().sum() == 24


def test_bootstrap_seed_reproduces_exposures(tmp_path):
    demo = tmp_path / "demo"
    assert run("demo", "--output-dir", str(demo)).returncode == 0
    outputs = []
    for index in range(2):
        output = tmp_path / f"bootstrap{index}.tsv"
        result = run(
            "signature",
            "-f",
            "MATRIX",
            "-i",
            str(demo / "matrix.tsv"),
            "-s",
            "5",
            "--bootstrap",
            "-br",
            "10",
            "--seed",
            "42",
            "-o",
            str(output),
        )
        assert result.returncode == 0, result.stderr
        outputs.append(output.read_text())
    assert outputs[0] == outputs[1]
