"""Check compressed input and sample identity through the profile command."""

import gzip
import subprocess
import sys

import pandas as pd

from tests.cli.local_fixtures import LOCAL_GENOME, LOCAL_MAF


def run(*args):
    return subprocess.run([sys.executable, "-m", "mutagene", *args], capture_output=True, text=True)


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
