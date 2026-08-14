"""Allele column resolution in MAF files (GitHub issue #98).

A MAF carrying Tumor_Seq_Allele2 without Tumor_Seq_Allele1 was read by the
profile parser but rejected by the context-window parser, so `profile` succeeded
on files where `rank`, `motif` and the webapp raised "Variant allele is not
defined in MAF file". These tests pin both parsers to the same acceptance rules.
"""

import io
from collections import namedtuple

import pytest

from mutagene.io import context_window, mutations_profile
from mutagene.io.maf_columns import resolve_alleles


def maf(**allele_columns):
    """A one-row MAF whose allele columns are exactly the ones given.

    Reference allele is C throughout; callers vary only the tumour columns.
    """
    header = ["Chromosome", "Start_Position", "End_Position", "Reference_Allele"]
    row = ["17", "7578406", "7578406", "C"]
    for name, value in allele_columns.items():
        header.append(name)
        row.append(value)
    header.append("Tumor_Sample_Barcode")
    row.append("SAMPLE1")
    return "\t".join(header) + "\n" + "\t".join(row) + "\n"


@pytest.fixture
def stub_genome(monkeypatch):
    """Give both parsers a fixed A/G context so no genome file is needed."""

    def window_context(mutations, twobit_file, window_size):
        return {
            (chrom, pos): (("A", "G"), [(chrom, pos, x, strand)])
            for chrom, pos, strand, x, y in mutations
        }

    def batch_context(mutations, assembly, method="twobit"):
        return {(chrom, pos): ("A", "G") for chrom, pos, _x, _y in mutations}

    monkeypatch.setattr(context_window, "get_context_twobit_window", window_context)
    monkeypatch.setattr(mutations_profile, "get_context_batch", batch_context)


def read_context_window(text):
    """Parse via the reader used by rank, motif and the webapp."""
    return context_window.read_MAF_with_context_window(io.StringIO(text), "unused.2bit", 1)


def read_profile(text):
    """Parse via the reader used by `mutagene profile`."""
    return mutations_profile.read_auto_profile(io.StringIO(text), fmt="MAF", asm=None)


class TestResolveAlleles:
    """Unit-level contract of the shared resolver."""

    def make(self, **columns):
        row = namedtuple("MAF", sorted(columns))
        return row(**columns)

    def test_allele2_only(self):
        data = self.make(reference_allele="C", tumor_seq_allele2="T")
        assert resolve_alleles(data) == ("C", "T")

    def test_allele1_only(self):
        data = self.make(reference_allele="C", tumor_seq_allele1="T")
        assert resolve_alleles(data) == ("C", "T")

    def test_allele2_wins_when_both_differ(self):
        data = self.make(reference_allele="C", tumor_seq_allele1="A", tumor_seq_allele2="T")
        assert resolve_alleles(data) == ("C", "T")

    def test_allele1_used_when_allele2_repeats_reference(self):
        data = self.make(reference_allele="C", tumor_seq_allele1="T", tumor_seq_allele2="C")
        assert resolve_alleles(data) == ("C", "T")

    def test_variant_allele_column_takes_precedence(self):
        data = self.make(reference_allele="C", variant_allele="G", tumor_seq_allele2="T")
        assert resolve_alleles(data) == ("C", "G")

    def test_no_variant_when_every_allele_matches_reference(self):
        data = self.make(reference_allele="C", tumor_seq_allele1="C", tumor_seq_allele2="C")
        assert resolve_alleles(data) == ("C", None)

    def test_missing_reference_raises(self):
        data = self.make(tumor_seq_allele2="T")
        with pytest.raises(ValueError, match="Reference allele"):
            resolve_alleles(data)

    def test_missing_every_variant_column_raises(self):
        data = self.make(reference_allele="C")
        with pytest.raises(ValueError, match="Variant allele"):
            resolve_alleles(data)


class TestContextWindowReader:
    """The reader that raised on the real 37-column file in #98."""

    def test_allele2_without_allele1(self, stub_genome):
        mutations, _, stats = read_context_window(maf(Tumor_Seq_Allele2="T"))
        assert stats["loaded"] == 1
        assert mutations["SAMPLE1"]["AGCT"] == 1.0

    def test_allele1_without_allele2(self, stub_genome):
        mutations, _, stats = read_context_window(maf(Tumor_Seq_Allele1="T"))
        assert stats["loaded"] == 1
        assert mutations["SAMPLE1"]["AGCT"] == 1.0

    def test_both_alleles_still_prefers_allele2(self, stub_genome):
        text = maf(Tumor_Seq_Allele1="C", Tumor_Seq_Allele2="T")
        mutations, _, stats = read_context_window(text)
        assert stats["loaded"] == 1
        assert mutations["SAMPLE1"]["AGCT"] == 1.0

    def test_no_allele_column_at_all_still_raises(self, stub_genome):
        with pytest.raises(ValueError, match="Variant allele is not defined"):
            read_context_window(maf())


class TestBothReadersAgree:
    """Neither parser may accept a header the other one rejects."""

    @pytest.mark.parametrize(
        "columns",
        [
            {"Tumor_Seq_Allele2": "T"},
            {"Tumor_Seq_Allele1": "T"},
            {"Tumor_Seq_Allele1": "C", "Tumor_Seq_Allele2": "T"},
            {"Variant_Allele": "T"},
        ],
        ids=["allele2-only", "allele1-only", "both-alleles", "variant-allele"],
    )
    def test_same_mutation_loaded_by_both(self, stub_genome, columns):
        text = maf(**columns)

        per_sample, _, window_stats = read_context_window(text)
        pooled, profile_stats = read_profile(text)

        assert window_stats["loaded"] == profile_stats["loaded"] == 1
        assert per_sample["SAMPLE1"]["AGCT"] == pooled["AGCT"] == 1.0
