"""Tests for tolerating malformed rows in MAF files (issues #92 and #61)."""

import logging
import os

import pytest

from mutagene.io import context_window
from mutagene.io.context_window import read_MAF_with_context_window

HEADER = (
    "Chromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele1\t"
    "Tumor_Seq_Allele2\tTumor_Sample_Barcode\tTranscript_Strand"
)

REAL_MAF = os.path.join(os.path.dirname(__file__), "..", "motifs", "data", "tcga_A3JM-01.maf")
REAL_GENOME = "test-reports/hg19.2bit"


def row(pos, strand="+", ref="C", alt="T", sample="SAMPLE1"):
    return f"17\t{pos}\t{ref}\t{ref}\t{alt}\t{sample}\t{strand}"


def maf(*rows):
    return [HEADER + "\n"] + [r + "\n" for r in rows]


@pytest.fixture
def fake_context(monkeypatch):
    """Replace the 2bit lookup so that parsing can be tested without a genome."""

    def get_context(mutations, twobit_file, window_size):
        return {
            (chrom, pos): (("A", "A"), [(chrom, pos, x, strand)])
            for chrom, pos, strand, x, y in mutations
        }

    monkeypatch.setattr(context_window, "get_context_twobit_window", get_context)


def read(lines):
    return read_MAF_with_context_window(lines, "unused.2bit", 1)


class TestWrongFieldCount:
    def test_keeps_valid_rows(self, fake_context):
        # Line 3 has two fields instead of seven
        lines = maf(row(100), "17\t200", row(300))
        mutations, mutations_with_context, stats = read(lines)

        assert stats["loaded"] == 2
        assert stats["skipped"] == 1
        assert [m[1] for m in mutations_with_context["SAMPLE1"]] == [100, 300]

    def test_warns_once_with_line_number(self, fake_context, caplog):
        lines = maf(row(100), "17\t200", row(300))
        with caplog.at_level(logging.WARNING, logger=context_window.__name__):
            read(lines)

        warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warnings) == 1
        assert "line 3" in warnings[0]
        assert "expected 7 fields, got 2" in warnings[0]

    def test_reports_every_skipped_row_in_a_single_warning(self, fake_context, caplog):
        lines = maf(row(100), "17\t200", row(300), "17", row(500))
        with caplog.at_level(logging.WARNING, logger=context_window.__name__):
            _, _, stats = read(lines)

        warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warnings) == 1
        assert "Skipped 2 malformed rows" in warnings[0]
        assert "line 3" in warnings[0] and "line 5" in warnings[0]
        assert stats["skipped"] == 2
        assert stats["loaded"] == 3

    def test_extra_fields_are_skipped_too(self, fake_context):
        lines = maf(row(100), row(200) + "\textra", row(300))
        _, _, stats = read(lines)

        assert stats["loaded"] == 2
        assert stats["skipped"] == 1

    def test_comment_lines_do_not_shift_reported_line_number(self, fake_context, caplog):
        lines = ["# a comment\n"] + maf(row(100), "17\t200")
        with caplog.at_level(logging.WARNING, logger=context_window.__name__):
            read(lines)

        assert "line 4" in caplog.records[0].message

    def test_blank_lines_are_not_counted_as_malformed(self, fake_context, caplog):
        lines = maf(row(100), "", row(300))
        with caplog.at_level(logging.WARNING, logger=context_window.__name__):
            _, _, stats = read(lines)

        assert stats["loaded"] == 2
        assert stats["skipped"] == 0
        assert caplog.records == []


class TestTranscriptStrand:
    def test_unexpected_value_skips_only_that_row(self, fake_context, caplog):
        lines = maf(row(100), row(200, strand="?"), row(300))
        with caplog.at_level(logging.WARNING, logger=context_window.__name__):
            _, mutations_with_context, stats = read(lines)

        assert stats["loaded"] == 2
        assert stats["skipped"] == 1
        assert [m[1] for m in mutations_with_context["SAMPLE1"]] == [100, 300]

        warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warnings) == 1
        assert "line 3" in warnings[0]
        assert "transcript_strand" in warnings[0]

    @pytest.mark.parametrize(
        "value,expected",
        [("+", "+"), ("-", "-"), ("1", "+"), ("-1", "-"), ("", "+")],
    )
    def test_known_values_are_normalized(self, fake_context, value, expected):
        _, mutations_with_context, stats = read(maf(row(100, strand=value)))

        assert stats["skipped"] == 0
        assert mutations_with_context["SAMPLE1"][0][2] == expected


class TestWellFormedInput:
    def test_no_warning_and_no_skips(self, fake_context, caplog):
        lines = maf(row(100), row(200, strand="-"), row(300, strand="1"))
        with caplog.at_level(logging.WARNING, logger=context_window.__name__):
            mutations, mutations_with_context, stats = read(lines)

        assert caplog.records == []
        assert stats == {"loaded": 3, "skipped": 0, "nsamples": 1, "format": "MAF"}
        assert sum(mutations["SAMPLE1"].values()) == 3

    def test_missing_transcript_strand_column_defaults_to_plus(self, fake_context):
        header = HEADER.rsplit("\t", 1)[0]
        lines = [header + "\n", row(100).rsplit("\t", 1)[0] + "\n"]
        _, mutations_with_context, stats = read(lines)

        assert stats["loaded"] == 1
        assert mutations_with_context["SAMPLE1"][0][2] == "+"


@pytest.mark.skipif(not os.path.isfile(REAL_GENOME), reason="hg19.2bit fixture not available")
def test_real_maf_with_short_row_completes():
    # Line 14 of this fixture has 83 fields instead of 95 (issue #92)
    with open(REAL_MAF) as f:
        _, mutations_with_context, stats = read_MAF_with_context_window(f, REAL_GENOME, 50)

    assert stats["loaded"] == 236
    assert stats["skipped"] == 1
    assert sum(len(v) for v in mutations_with_context.values()) == 236
