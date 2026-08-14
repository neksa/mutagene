"""Flanking-base handling when REF is on the reverse strand (GitHub issue #101).

`get_context_53_twobit` complemented the flanking bases one at a time, so the
second assignment read the value the first had just written and the 5' base was
never complemented at all. `get_context_twobit_window` does the same swap
correctly, so the two context extractors disagreed for the same input.

These tests pin the two against each other rather than only against literal
expected values, so the pair cannot drift apart again.
"""

import pytest
import twobitreader

from mutagene.dna import complementary_nucleotide as cn
from mutagene.io.context_window import get_context_twobit_window
from mutagene.io.mutations_profile import get_context_53_twobit


class FakeChromosome:
    def __init__(self, sequence):
        self.sequence = sequence

    def __getitem__(self, key):
        return self.sequence[key]


class FakeTwoBitFile:
    """Minimal stand-in for twobitreader.TwoBitFile over in-memory sequences."""

    def __init__(self, sequences):
        self.sequences = sequences

    def __contains__(self, name):
        return name in self.sequences

    def __getitem__(self, name):
        return FakeChromosome(self.sequences[name])


@pytest.fixture
def genome(monkeypatch):
    """A chr17 whose base 2 is G, flanked by A on both sides.

    Position 2 is the mutated base. A MAF reporting REF=C there is describing
    the reverse strand, which is the branch under test. Both flanks are A so
    that the correct answer (T, T) differs from what the buggy swap produced.
    """
    fake = FakeTwoBitFile({"chr17": "AGA"})
    monkeypatch.setattr(twobitreader, "TwoBitFile", lambda fname: fake)
    return fake


POSITION = 2
REVERSE_STRAND_REF = "C"  # complement of the G actually in the assembly


def test_reverse_strand_context_is_complemented_and_swapped(genome):
    contexts = get_context_53_twobit([("17", POSITION, REVERSE_STRAND_REF, "T")], "unused")

    assert contexts[("17", POSITION)] == (cn["A"], cn["A"]) == ("T", "T")


def test_both_context_extractors_agree_on_reverse_strand(genome):
    """The profile reader and the context-window reader must not diverge."""
    mutation = ("17", POSITION, REVERSE_STRAND_REF, "T")

    from_profile = get_context_53_twobit([mutation], "unused")[("17", POSITION)]

    windowed = get_context_twobit_window(
        [("17", POSITION, "+", REVERSE_STRAND_REF, "T")], "unused", 1
    )
    (nuc5, nuc3), _seq = windowed[("17", POSITION)]

    assert from_profile == (nuc5, nuc3)


def test_forward_strand_context_is_untouched(genome):
    """When REF matches the assembly the flanks are returned as they are."""
    contexts = get_context_53_twobit([("17", POSITION, "G", "A")], "unused")

    assert contexts[("17", POSITION)] == ("A", "A")


def test_ref_matching_neither_strand_yields_no_context(genome):
    contexts = get_context_53_twobit([("17", POSITION, "T", "A")], "unused")

    assert contexts[("17", POSITION)] == ("N", "N")
