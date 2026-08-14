"""One unusable sample must not discard the rest of the file (issue #100).

The per-sample loop returned from the whole function when `get_context_twobit_window`
came back empty, so a single sample whose chromosomes were missing from the
assembly threw away every sample after it and the caller saw `loaded: 0` as
though the file were empty. Sample order decided how much data was lost.
"""

import logging

import pytest

from mutagene.io import context_window

HEADER = "Chromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode"


def row(chrom, sample, pos=100):
    return f"{chrom}\t{pos}\tC\tT\t{sample}"


def maf(*rows):
    return [HEADER + "\n"] + [r + "\n" for r in rows]


@pytest.fixture
def genome_with_only_chr17(monkeypatch):
    """An assembly that knows chr17, so mutations elsewhere yield no context."""

    def get_context(mutations, twobit_file, window_size):
        return {
            (chrom, pos): (("A", "G"), [(chrom, pos, x, strand)])
            for chrom, pos, strand, x, y in mutations
            if chrom == "17"
        }

    monkeypatch.setattr(context_window, "get_context_twobit_window", get_context)


def read(lines):
    return context_window.read_MAF_with_context_window(lines, "unused.2bit", 1)


class TestSampleWithoutContext:
    def test_later_samples_survive_an_unusable_first_sample(self, genome_with_only_chr17):
        lines = maf(row("99", "BAD_SAMPLE"), row("17", "GOOD_SAMPLE"))
        mutations, mutations_with_context, stats = read(lines)

        assert stats["loaded"] == 1, "the good sample was discarded with the bad one"
        assert stats["nsamples"] == 1
        assert "GOOD_SAMPLE" in mutations
        assert "BAD_SAMPLE" not in mutations
        assert [m[1] for m in mutations_with_context["GOOD_SAMPLE"]] == [100]

    def test_result_does_not_depend_on_sample_order(self, genome_with_only_chr17):
        good_first = read(maf(row("17", "GOOD_SAMPLE"), row("99", "BAD_SAMPLE")))
        bad_first = read(maf(row("99", "BAD_SAMPLE"), row("17", "GOOD_SAMPLE")))

        assert good_first[0] == bad_first[0]
        assert good_first[2]["loaded"] == bad_first[2]["loaded"] == 1

    def test_unusable_mutations_are_counted_as_skipped(self, genome_with_only_chr17):
        lines = maf(row("99", "BAD_SAMPLE"), row("99", "BAD_SAMPLE", pos=200), row("17", "OK"))
        _, _, stats = read(lines)

        assert stats["skipped"] == 2
        assert stats["loaded"] == 1

    def test_warns_once_naming_the_dropped_samples(self, genome_with_only_chr17, caplog):
        lines = maf(row("99", "BAD_SAMPLE"), row("17", "GOOD_SAMPLE"))
        with caplog.at_level(logging.WARNING, logger=context_window.__name__):
            read(lines)

        warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warnings) == 1
        assert "BAD_SAMPLE" in warnings[0]
        assert "GOOD_SAMPLE" not in warnings[0]

    def test_every_sample_unusable_reports_an_empty_result_not_none(self, genome_with_only_chr17):
        """Callers index into these, so returning None turned a bad file into a crash."""
        mutations, mutations_with_context, stats = read(maf(row("99", "BAD_SAMPLE")))

        assert mutations == {}
        assert mutations_with_context == {}
        assert stats["loaded"] == 0
        assert stats["format"] == "MAF"

    def test_no_warning_when_every_sample_resolves(self, genome_with_only_chr17, caplog):
        lines = maf(row("17", "A"), row("17", "B"))
        with caplog.at_level(logging.WARNING, logger=context_window.__name__):
            _, _, stats = read(lines)

        assert stats["loaded"] == 2
        assert stats["nsamples"] == 2
        assert caplog.records == []


class TestOtherReaders:
    """The same abort was present in the TCGI and VCF readers."""

    def test_tcgi_keeps_usable_samples(self, genome_with_only_chr17):
        lines = [
            "CHROM\tPOS\tREF\tALT\tSAMPLE\n",
            "99\t100\tC\tT\tBAD_SAMPLE\n",
            "17\t100\tC\tT\tGOOD_SAMPLE\n",
        ]
        mutations, _, stats = context_window.read_TCGI_with_context_window(lines, "unused.2bit", 1)

        assert stats["loaded"] == 1
        assert list(mutations) == ["GOOD_SAMPLE"]

    def test_vcf_without_context_returns_empty_not_none(self, genome_with_only_chr17):
        lines = [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "chr99\t100\t.\tC\tT\t.\tPASS\t.\n",
        ]
        mutations, mutations_with_context, stats = context_window.read_VCF_with_context_window(
            lines, "unused.2bit", 1
        )

        assert mutations == {}
        assert mutations_with_context == {}
        assert stats["loaded"] == 0
