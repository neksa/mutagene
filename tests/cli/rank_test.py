import json

import pandas as pd

from tests.cli import cli_test_utils

DEFAULT_INFILE = f"{cli_test_utils.TEST_DIR}/sample1.maf"
DEFAULT_GENOME = f"{cli_test_utils.TEST_DIR}/hg19.2bit"


def read_ranking(outfile):
    """Read a ranking the way the documentation says to: skipping the provenance."""
    return pd.read_csv(outfile, sep="\t", comment="#")


def read_provenance_header(outfile):
    with open(outfile) as fh:
        return [line[2:].strip() for line in fh if line.startswith("# ")]


# mutagene -v rank -g hg19 -i sample1.maf -c pancancer -o test-reports/cli-rank-sample1-pancancer.txt
def test_rank(test_data):
    infile = DEFAULT_INFILE
    outfile = f"{cli_test_utils.TEST_DIR}/cli-rank-sample1-pancancer.txt"
    genome = DEFAULT_GENOME

    cli_test_utils.run_with_args(
        "rank", ["-i", infile, "-o", outfile, "-g", genome, "-c", "pancancer"]
    )

    ranking = read_ranking(outfile)

    assert ranking.iloc[0]["gene"] == "CPXM2"
    assert ranking.iloc[0]["mutation"] == "T536M"


# mutagene rank -i sample1.maf -g hg19 -c gcb_lymphomas -o test-reports/cli-rank-sample1-gcb_lymphomas.txt
def test_rank_4_1(test_data):
    infile = DEFAULT_INFILE
    outfile = f"{cli_test_utils.TEST_DIR}/cli-rank-sample1-gcb_lymphomas.txt"
    genome = DEFAULT_GENOME

    cli_test_utils.run_with_args(
        "rank", ["-i", infile, "-o", outfile, "-g", genome, "-c", "gcb_lymphomas"]
    )

    ranking = read_ranking(outfile)

    assert ranking.iloc[0]["gene"] == "BOD1L"
    assert ranking.iloc[0]["mutation"] == "T2810S"


class TestProvenance:
    """The output must say what produced it (GitHub issue #63)."""

    def outfile(self):
        return f"{cli_test_utils.TEST_DIR}/cli-rank-provenance.txt"

    def run(self, extra_args=()):
        outfile = self.outfile()
        cli_test_utils.run_with_args(
            "rank",
            ["-i", DEFAULT_INFILE, "-o", outfile, "-g", DEFAULT_GENOME, "-c", "pancancer"]
            + list(extra_args),
        )
        return outfile

    def test_header_names_every_input_source(self, test_data):
        header = "\n".join(read_provenance_header(self.run()))

        assert "profile_source: precalculated cohort pancancer" in header
        assert "observed_mutations_source: precalculated cohort pancancer" in header
        assert "cohort_size_source: precalculated cohort pancancer" in header
        assert "cohort_size:" in header
        assert "threshold_driver:" in header and "threshold_passenger:" in header
        assert "mutagene_version:" in header

    def test_sidecar_matches_the_header(self, test_data):
        outfile = self.run()

        with open(f"{outfile}.provenance.json") as fh:
            sidecar = json.load(fh)

        header = read_provenance_header(outfile)
        assert [f"{k}: {v}" for k, v in sidecar.items()] == header

    def test_nsamples_override_is_recorded_as_the_source(self, test_data):
        header = "\n".join(read_provenance_header(self.run(["-n", "42"])))

        assert "cohort_size: 42" in header
        assert "cohort_size_source: --nsamples" in header

    def test_table_is_still_machine_readable(self, test_data):
        ranking = read_ranking(self.run())

        assert list(ranking.columns) == [
            "gene",
            "transcript",
            "mutation",
            "mutability",
            "observed",
            "bscore",
            "qvalue",
            "label",
        ]
        assert len(ranking) > 0
        assert ranking["bscore"].is_monotonic_increasing


class TestTranscriptSelection:
    """One transcript per gene, chosen the same way every run (GitHub issue #48)."""

    def test_every_gene_uses_a_single_transcript(self, test_data):
        outfile = f"{cli_test_utils.TEST_DIR}/cli-rank-transcripts.txt"
        cli_test_utils.run_with_args(
            "rank", ["-i", DEFAULT_INFILE, "-o", outfile, "-g", DEFAULT_GENOME, "-c", "pancancer"]
        )

        ranking = read_ranking(outfile)
        per_gene = ranking.groupby("gene")["transcript"].nunique()

        assert per_gene.max() == 1, "a gene was ranked under more than one transcript"
        assert ranking["transcript"].notna().all()
