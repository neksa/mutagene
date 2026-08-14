"""Detecting a wrong genome assembly (GitHub issue #99).

The webapp used to count log records emitted by `context_window`, attached only
for the duration of `calc_profile`. That could never work: the message it looked
for comes from a different code path and is emitted after the handler is
removed, and the path `calc_profile` actually takes never logged it at all. So
analyses run against the wrong assembly completed silently.

The counters now come from the code that compares the alleles.
"""

import io

import pytest
import twobitreader

from mutagene.io import mutations_profile
from mutagene.io.context_stats import (
    MISMATCH_MIN_COUNT,
    assembly_mismatch_warning,
    merge_context_stats,
    mismatch_rate,
    new_context_stats,
)
from mutagene.io.context_window import get_context_twobit_window, read_MAF_with_context_window
from mutagene.io.mutations_profile import get_context_53_twobit

HEADER = (
    "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
    "Tumor_Seq_Allele2\tTumor_Sample_Barcode"
)


class FakeChromosome:
    def __init__(self, sequence):
        self.sequence = sequence

    def __getitem__(self, key):
        return self.sequence[key]


class FakeTwoBitFile:
    def __init__(self, sequences):
        self.sequences = sequences

    def __contains__(self, name):
        return name in self.sequences

    def __getitem__(self, name):
        return FakeChromosome(self.sequences[name])


@pytest.fixture
def assembly(monkeypatch):
    """A chr17 that is all A, so any non-A/T reference allele is a mismatch."""
    fake = FakeTwoBitFile({"chr17": "A" * 200})
    monkeypatch.setattr(twobitreader, "TwoBitFile", lambda fname: fake)
    return fake


def maf(count, ref="G"):
    """A MAF whose reference allele is `ref` at consecutive positions."""
    rows = [HEADER]
    for i in range(count):
        pos = 10 + i
        rows.append(f"17\t{pos}\t{pos}\t{ref}\tT\tSAMPLE1")
    return [line + "\n" for line in rows]


class TestCountingAtTheSource:
    def test_window_extractor_counts_every_mismatch(self, assembly):
        # REF=G against an all-A assembly matches neither strand.
        mutations = [("17", 10 + i, "+", "G", "T") for i in range(7)]

        _contexts, stats = get_context_twobit_window(mutations, "unused", 1)

        assert stats["ref_mismatches"] == 7
        assert stats["reverse_strand_ref"] == 0

    def test_profile_extractor_counts_them_too(self, assembly):
        """This is the path calc_profile takes; it reported nothing at all before."""
        mutations = [("17", 10 + i, "G", "C") for i in range(7)]

        _contexts, stats = get_context_53_twobit(mutations, "unused")

        assert stats["ref_mismatches"] == 7

    def test_reverse_strand_ref_is_not_an_assembly_mismatch(self, assembly):
        """REF=T against an A is the complementary strand, which is handled, not wrong."""
        mutations = [("17", 10, "+", "T", "C")]

        _contexts, stats = get_context_twobit_window(mutations, "unused", 1)

        assert stats["ref_mismatches"] == 0
        assert stats["reverse_strand_ref"] == 1

    def test_matching_reference_produces_no_counts(self, assembly):
        mutations = [("17", 10, "+", "A", "C")]

        _contexts, stats = get_context_twobit_window(mutations, "unused", 1)

        assert stats == new_context_stats()

    def test_unknown_chromosome_is_counted_separately(self, assembly):
        mutations = [("99", 10, "+", "A", "C")]

        _contexts, stats = get_context_twobit_window(mutations, "unused", 1)

        assert stats["chromosome_not_found"] == 1
        assert stats["ref_mismatches"] == 0


class TestStatsReachTheCaller:
    def test_read_mutations_reports_mismatches(self, assembly):
        _mutations, _ctx, stats = read_MAF_with_context_window(maf(5), "unused", 1)

        assert stats["ref_mismatches"] == 5
        assert stats["loaded"] == 0, "mismatched mutations have no usable context"

    def test_read_profile_reports_mismatches(self, assembly):
        text = "".join(maf(5))
        _mutations, stats = mutations_profile.read_auto_profile(
            io.StringIO(text), fmt="MAF", asm="unused"
        )

        assert stats["ref_mismatches"] == 5

    def test_matching_assembly_reports_none(self, assembly):
        _mutations, _ctx, stats = read_MAF_with_context_window(maf(5, ref="A"), "unused", 1)

        assert stats["ref_mismatches"] == 0
        assert stats["loaded"] == 5


class TestWarningThreshold:
    def rate_stats(self, mismatches):
        stats = new_context_stats()
        stats["ref_mismatches"] = mismatches
        return stats

    def test_no_warning_below_the_count_threshold(self):
        assert assembly_mismatch_warning(self.rate_stats(MISMATCH_MIN_COUNT), 0, "hg19") is None

    def test_warning_when_most_mutations_mismatch(self):
        warning = assembly_mismatch_warning(self.rate_stats(500), 100, "hg19")

        assert warning is not None
        assert warning["mismatch_count"] == 500
        assert "hg38" in warning["message"], "should name the other assembly"

    def test_no_warning_when_mismatches_are_a_rounding_error(self):
        """A handful of odd rows in a huge file is not a wrong assembly."""
        assert assembly_mismatch_warning(self.rate_stats(25), 4_000_000, "hg19") is None

    def test_rate_is_over_everything_examined(self):
        assert mismatch_rate(self.rate_stats(50), 50) == pytest.approx(0.5)

    def test_rate_of_a_file_where_nothing_survived(self):
        assert mismatch_rate(self.rate_stats(10), 0) == pytest.approx(1.0)

    def test_rate_without_any_mutations_is_zero_not_an_error(self):
        assert mismatch_rate(new_context_stats(), 0) == 0.0

    def test_unknown_assembly_name_omits_the_suggestion(self):
        warning = assembly_mismatch_warning(self.rate_stats(500), 100, "/genomes/custom.2bit")

        assert warning is not None
        assert "Please try" not in warning["message"]
        assert warning["genome"] == "custom"


class TestMergeContextStats:
    def test_sums_each_counter(self):
        total = new_context_stats()
        total["ref_mismatches"] = 2
        merge_context_stats(total, {"ref_mismatches": 3, "read_errors": 1})

        assert total["ref_mismatches"] == 5
        assert total["read_errors"] == 1

    def test_ignores_unrelated_keys(self):
        total = new_context_stats()
        merge_context_stats(total, {"loaded": 99, "format": "MAF"})

        assert total == new_context_stats()
