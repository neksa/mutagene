"""The analysis commands, exercised without a network (GitHub issue #125).

The CLI tests used a genome downloaded from UCSC, so a runner that could not
reach it ran fewer tests and reported ~5% lower coverage for the same commit.
These use a committed genome and a MAF generated to match it.
"""

import pandas as pd
import pytest

from tests.cli import cli_test_utils
from tests.cli.local_fixtures import (
    LOCAL_GENOME,
    LOCAL_MAF,
    LOCAL_MAF_MUTATIONS,
    LOCAL_MAF_SAMPLES,
)


def run(command, args):
    cli_test_utils.run_with_args(command, args)


class TestProfile:
    def test_writes_ninety_six_channels(self, tmp_path):
        out = tmp_path / "profile.tsv"
        run("profile", ["-i", LOCAL_MAF, "-g", LOCAL_GENOME, "-o", str(out)])

        rows = [line for line in out.read_text().splitlines() if line.strip()]
        assert len(rows) == 96
        assert sum(int(r.split("\t")[1]) for r in rows) == LOCAL_MAF_MUTATIONS

    def test_every_reference_allele_matches_the_fixture_genome(self, tmp_path, caplog):
        import logging

        out = tmp_path / "profile.tsv"
        with caplog.at_level(logging.WARNING):
            run("profile", ["-i", LOCAL_MAF, "-g", LOCAL_GENOME, "-o", str(out)])

        assert not [r for r in caplog.records if "matches neither strand" in r.message]


class TestSignature:
    @pytest.mark.parametrize("collection", ["MGA", "COSMICv2"])
    def test_decomposes(self, tmp_path, collection):
        out = tmp_path / "sig.tsv"
        run(
            "signature",
            ["-i", LOCAL_MAF, "-g", LOCAL_GENOME, "-s", collection, "-o", str(out)],
        )

        assert out.stat().st_size > 0


class TestMotif:
    def test_a_custom_motif_runs(self, tmp_path):
        out = tmp_path / "motif.tsv"
        run("motif", ["-i", LOCAL_MAF, "-g", LOCAL_GENOME, "--motif", "C[C>T]G", "-o", str(out)])

        assert out.exists()

    def test_the_builtin_motifs_run(self, tmp_path):
        out = tmp_path / "motif.tsv"
        run("motif", ["-i", LOCAL_MAF, "-g", LOCAL_GENOME, "-o", str(out)])

        assert out.exists()


class TestRankWithoutCohort:
    def test_ranking_the_input_against_itself_reports_no_annotations(self, tmp_path, caplog):
        """The fixture has no protein annotations, which rank should say plainly."""
        import logging

        out = tmp_path / "rank.tsv"
        with caplog.at_level(logging.WARNING):
            run("rank", ["-i", LOCAL_MAF, "-g", LOCAL_GENOME, "-o", str(out)])

        assert any("could be ranked" in r.message for r in caplog.records)


class TestParamsRoundTrip:
    def test_a_recorded_run_reproduces_its_output(self, tmp_path):
        first, second = tmp_path / "a.tsv", tmp_path / "b.tsv"
        params = tmp_path / "run.json"

        run(
            "profile",
            ["-i", LOCAL_MAF, "-g", LOCAL_GENOME, "-o", str(first), "--params-out", str(params)],
        )
        run("profile", ["--params-in", str(params), "-o", str(second)])

        assert first.read_text() == second.read_text()


class TestFilterEndToEnd:
    def test_filtered_variants_are_excluded(self, tmp_path):
        """sample_small.maf carries a FILTER column, all PASS."""
        rows = LOCAL_MAF and open(LOCAL_MAF).read().splitlines()
        rejected = tmp_path / "rejected.maf"
        rejected.write_text(
            "\n".join([rows[0]] + [r.replace("\tPASS", "\tgermline") for r in rows[1:]]) + "\n"
        )

        kept = tmp_path / "kept.tsv"
        with pytest.raises(SystemExit):
            run("profile", ["-i", str(rejected), "-g", LOCAL_GENOME, "-o", str(kept)])

    def test_keep_filtered_reads_them_again(self, tmp_path):
        rows = open(LOCAL_MAF).read().splitlines()
        rejected = tmp_path / "rejected.maf"
        rejected.write_text(
            "\n".join([rows[0]] + [r.replace("\tPASS", "\tgermline") for r in rows[1:]]) + "\n"
        )

        out = tmp_path / "kept.tsv"
        run("profile", ["-i", str(rejected), "-g", LOCAL_GENOME, "-o", str(out), "--keep-filtered"])

        total = sum(
            int(line.split("\t")[1]) for line in out.read_text().splitlines() if line.strip()
        )
        assert total == LOCAL_MAF_MUTATIONS
