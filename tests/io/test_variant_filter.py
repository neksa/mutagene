"""Honouring the caller's FILTER field.

FILTER is a VCF column that GDC-style MAFs carry forward, and neither reader
looked at it. Reading rejected variants as accepted changes what an analysis is
about: a Mutect2 file where 4,637 of 4,675 rows are marked `germline` gives a
"mutational profile" that is mostly germline variation, and the decomposition
duly reports germline-contamination signatures.
"""

import io

import pytest

from mutagene.io import context_window, mutations_profile
from mutagene.io.context_stats import new_context_stats
from mutagene.io.variant_filter import passes_filter

MAF_HEADER = (
    "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
    "Tumor_Seq_Allele2\tTumor_Sample_Barcode\tFILTER"
)
VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"


@pytest.fixture
def stub_genome(monkeypatch):
    def window(mutations, twobit_file, window_size):
        contexts = {
            (chrom, pos): (("A", "G"), [(chrom, pos, x, strand)])
            for chrom, pos, strand, x, y in mutations
        }
        return contexts, new_context_stats()

    def batch(mutations, assembly, method="twobit"):
        return {(chrom, pos): ("A", "G") for chrom, pos, _x, _y in mutations}, new_context_stats()

    monkeypatch.setattr(context_window, "get_context_twobit_window", window)
    monkeypatch.setattr(mutations_profile, "get_context_batch", batch)


def maf(*filters):
    rows = [MAF_HEADER]
    for i, value in enumerate(filters):
        pos = 100 + i
        rows.append(f"17\t{pos}\t{pos}\tC\tT\tSAMPLE1\t{value}")
    return "\n".join(rows) + "\n"


def vcf(*filters):
    rows = [VCF_HEADER]
    for i, value in enumerate(filters):
        rows.append(f"chr17\t{100 + i}\t.\tC\tT\t.\t{value}\t.")
    return "\n".join(rows) + "\n"


class TestPassesFilter:
    @pytest.mark.parametrize("value", ["PASS", ".", "", "  PASS  "])
    def test_accepted_values(self, value):
        assert passes_filter(value) is True

    @pytest.mark.parametrize(
        "value", ["germline", "weak_evidence", "clustered_events;germline", "map_qual"]
    )
    def test_rejected_values(self, value):
        assert passes_filter(value) is False

    def test_a_missing_column_passes(self):
        """A file that records no filters has not rejected anything."""
        assert passes_filter(None) is True


class TestMafReaders:
    def test_profile_reader_keeps_only_passing_variants(self, stub_genome):
        _m, stats = mutations_profile.read_auto_profile(
            io.StringIO(maf("PASS", "germline", "PASS", "weak_evidence")), fmt="MAF", asm=None
        )

        assert stats["loaded"] == 2

    def test_context_window_reader_keeps_only_passing_variants(self, stub_genome):
        _m, _c, stats = context_window.read_MAF_with_context_window(
            io.StringIO(maf("PASS", "germline", "PASS", "weak_evidence")), "unused", 1
        )

        assert stats["loaded"] == 2

    def test_both_readers_agree(self, stub_genome):
        text = maf("PASS", "germline", ".", "clustered_events;germline")

        _m, profile_stats = mutations_profile.read_auto_profile(
            io.StringIO(text), fmt="MAF", asm=None
        )
        _m2, _c, window_stats = context_window.read_MAF_with_context_window(
            io.StringIO(text), "unused", 1
        )

        assert profile_stats["loaded"] == window_stats["loaded"] == 2

    def test_keep_filtered_includes_everything(self, stub_genome):
        _m, stats = mutations_profile.read_auto_profile(
            io.StringIO(maf("PASS", "germline", "weak_evidence")),
            fmt="MAF",
            asm=None,
            keep_filtered=True,
        )

        assert stats["loaded"] == 3

    def test_a_file_without_the_column_is_untouched(self, stub_genome):
        """Most MAFs have no FILTER column; they must not lose a single row."""
        text = (
            "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
            "Tumor_Seq_Allele2\tTumor_Sample_Barcode\n"
            "17\t100\t100\tC\tT\tSAMPLE1\n17\t200\t200\tC\tT\tSAMPLE1\n"
        )
        _m, stats = mutations_profile.read_auto_profile(io.StringIO(text), fmt="MAF", asm=None)

        assert stats["loaded"] == 2

    def test_gdc_filter_is_not_used(self, stub_genome):
        """GDC_FILTER mixes quality flags with annotations such as NonExonic,
        and a non-exonic mutation is perfectly good input to a profile."""
        text = (
            "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
            "Tumor_Seq_Allele2\tTumor_Sample_Barcode\tGDC_FILTER\n"
            "17\t100\t100\tC\tT\tSAMPLE1\tNonExonic\n"
        )
        _m, stats = mutations_profile.read_auto_profile(io.StringIO(text), fmt="MAF", asm=None)

        assert stats["loaded"] == 1


class TestVcfReaders:
    def test_profile_reader_keeps_only_passing_variants(self, stub_genome):
        _m, stats = mutations_profile.read_auto_profile(
            io.StringIO(vcf("PASS", "germline", ".", "map_qual")), fmt="VCF", asm=None
        )

        assert stats["loaded"] == 2

    def test_context_window_reader_keeps_only_passing_variants(self, stub_genome):
        lines = vcf("PASS", "germline", ".", "map_qual").splitlines(keepends=True)
        _m, _c, stats = context_window.read_VCF_with_context_window(lines, "unused", 1)

        assert stats["loaded"] == 2

    def test_keep_filtered_includes_everything(self, stub_genome):
        _m, stats = mutations_profile.read_auto_profile(
            io.StringIO(vcf("PASS", "germline", "map_qual")),
            fmt="VCF",
            asm=None,
            keep_filtered=True,
        )

        assert stats["loaded"] == 3

    def test_a_row_that_stops_before_filter_is_kept(self, stub_genome):
        """Five columns is a valid minimal VCF row; it records no filters."""
        text = "#CHROM\tPOS\tID\tREF\tALT\nchr17\t100\t.\tC\tT\n"
        _m, stats = mutations_profile.read_auto_profile(io.StringIO(text), fmt="VCF", asm=None)

        assert stats["loaded"] == 1
