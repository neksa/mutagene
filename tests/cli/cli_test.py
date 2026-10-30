import json
import os

import pytest

from tests.cli import cli_test_utils


def test_fetch():
    # mutagene -v fetch cohorts MSKCC --cohort paac_jhu_2014
    cli_test_utils.run_with_args("fetch", ["cohorts", "MSKCC", "--cohort", "paac_jhu_2014"])

    file_name = "paac_jhu_2014.tar.gz"
    os.rename(f"./{file_name}", f"{cli_test_utils.TEST_DIR}/{file_name}")

    # File size check: external data source may change, so check for reasonable size (> 100KB)
    file_size = os.path.getsize(f"{cli_test_utils.TEST_DIR}/{file_name}")
    assert file_size > 100000, f"Downloaded file too small: {file_size} bytes"


def test_profile(test_data):
    infile = "tests/motifs/data/vcf/data.vcf"
    outfile = f"{cli_test_utils.TEST_DIR}/cli-profile-test.txt"
    genome = f"{cli_test_utils.TEST_DIR}/hg19.2bit"  # 'tests/motifs/data/test_genome.2bit'

    # mutagene -v profile -i sample1.maf -g hg19 -o test-reports/cli-profile-sample1.txt
    cli_test_utils.run_with_args("profile", ["-i", infile, "-o", outfile, "-g", genome])

    out_lines = cli_test_utils.get_file_lines(outfile, 3)

    assert out_lines[0] == "A[C>A]A\t0\n"
    assert out_lines[1] == "A[C>G]A\t2\n"
    assert out_lines[2] == "A[C>T]A\t4\n"


def test_signature(test_data):
    infile = f"{cli_test_utils.TEST_DIR}/sample1.maf"
    outfile = f"{cli_test_utils.TEST_DIR}/cli-signature-sample1.txt"
    diagnostics = f"{cli_test_utils.TEST_DIR}/cli-signature-sample1.json"
    genome = f"{cli_test_utils.TEST_DIR}/hg19.2bit"

    # mutagene -v signature identify -i sample1.maf -g hg19 -s5 -o test-reports/cli-signature-sample1.txt
    cli_test_utils.run_with_args(
        "signature",
        [
            "identify",
            "-i",
            infile,
            "-o",
            outfile,
            "-g",
            genome,
            "-s5",
            "--diagnostics-out",
            diagnostics,
        ],
    )

    out_lines = cli_test_utils.get_file_lines(outfile, 2)

    fields = out_lines[1].strip().split("\t")
    assert fields[:2] == ["TCGA-50-6593-01A-11D-1753-08", "MGA-1"]
    # Sparse BIC selects components before the final constrained likelihood refit.
    assert float(fields[2]) == pytest.approx(0.544836, abs=1e-5)
    assert int(fields[3]) == 185
    with open(diagnostics) as stream:
        fit = json.load(stream)["samples"][fields[0]]
    assert fit["maximum_refit_optimality_gap"] <= 2e-5
    assert fit["background_selection"] == "joint"
    assert fit["estimator"] == "sparse_multinomial_mle"
    assert fit["fallback"] is None
