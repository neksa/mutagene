"""Unreadable input must produce a message, not a traceback.

Every case here raised an unhandled exception out of the parser: an IndexError
from a short VCF row, a ValueError from non-numeric text sniffed as VCF, an
AttributeError from a missing MAF column, and a TypeError in the caller because
an unrecognised file returned None where a stats dict was expected.

A bad *row* is skipped and counted. A bad *file* raises ValueError, which the
CLI turns into a message.
"""

import io

import pytest

from mutagene.io import mutations_profile
from mutagene.io.context_stats import new_context_stats


@pytest.fixture
def stub_context(monkeypatch):
    def fake(mutations, assembly, method="twobit"):
        contexts = {(chrom, pos): ("A", "G") for chrom, pos, _x, _y in mutations}
        return contexts, new_context_stats()

    monkeypatch.setattr(mutations_profile, "get_context_batch", fake)


def read(text, fmt=None):
    return mutations_profile.read_auto_profile(io.StringIO(text), fmt=fmt, asm=None)


class TestUnrecognisableFile:
    """These returned (None, None), which the caller then indexed into."""

    def test_empty_file_reports_empty_stats(self):
        mutations, stats = read("")

        assert mutations == {}
        assert stats["loaded"] == 0
        assert stats["format"] == "unknown"

    def test_prose_reports_empty_stats(self):
        mutations, stats = read("this is not a maf at all\njust prose\n")

        assert stats["loaded"] == 0

    def test_header_without_rows(self):
        mutations, stats = read("Chromosome\tStart_Position\n")

        assert stats["loaded"] == 0

    def test_stats_are_always_subscriptable(self):
        """calc_profile reads stats['loaded'] directly, so None was a TypeError."""
        for text in ("", "nonsense\n", "Chromosome\tStart_Position\n"):
            _mutations, stats = read(text)
            assert stats["loaded"] == 0 and stats["skipped"] == 0


class TestMalformedVcfRows:
    def test_row_with_no_alt_column_is_skipped(self, stub_context):
        """ALT is the fifth column; the guard asked for four and indexed the fifth."""
        _mutations, stats = read("chr17\t100\t.\tC\n", fmt="VCF")

        assert stats["skipped"] == 1
        assert stats["loaded"] == 0

    def test_non_numeric_position_is_skipped(self, stub_context):
        _mutations, stats = read("chr17\tnot_a_number\t.\tC\tT\n", fmt="VCF")

        assert stats["skipped"] == 1

    def test_good_rows_survive_bad_ones(self, stub_context):
        text = "chr17\t100\t.\tC\tT\nchr17\tbad\t.\tC\tT\nchr17\t200\t.\tC\tT\n"
        mutations, stats = read(text, fmt="VCF")

        assert stats["loaded"] == 2, "a bad row cost the good ones"
        assert stats["skipped"] == 1
        assert mutations["AGCT"] == 2.0

    def test_skipped_rows_are_reported_once(self, stub_context, caplog):
        import logging

        text = "chr17\tbad\t.\tC\tT\nchr17\talso_bad\t.\tC\tT\n"
        with caplog.at_level(logging.WARNING):
            read(text, fmt="VCF")

        warnings = [r.message for r in caplog.records if "malformed" in r.message]
        assert len(warnings) == 1
        assert "Skipped 2 malformed rows" in warnings[0]


class TestMalformedMafRows:
    HEADER = "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2"

    def maf(self, *rows):
        return "\n".join([self.HEADER, *rows]) + "\n"

    def test_non_numeric_position_is_skipped(self, stub_context):
        text = self.maf("17\tbad\tbad\tC\tT", "17\t100\t100\tC\tT")
        mutations, stats = read(text, fmt="MAF")

        assert stats["loaded"] == 1
        assert stats["skipped"] == 1

    def test_wrong_field_count_is_skipped(self, stub_context):
        text = self.maf("17\t100", "17\t200\t200\tC\tT")
        _mutations, stats = read(text, fmt="MAF")

        assert stats["loaded"] == 1
        assert stats["skipped"] == 1

    def test_missing_position_column_is_a_file_level_error(self, stub_context):
        """A missing column is a property of the file, so it stops the read."""
        text = "Chromosome\tReference_Allele\tTumor_Seq_Allele2\n17\tC\tT\n"

        with pytest.raises(ValueError, match="Start_Position and End_Position"):
            read(text, fmt="MAF")

    def test_missing_chromosome_column_is_a_file_level_error(self, stub_context):
        text = "Start_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\n1\t1\tC\tT\n"

        with pytest.raises(ValueError, match="Chromosome is not defined"):
            read(text, fmt="MAF")
