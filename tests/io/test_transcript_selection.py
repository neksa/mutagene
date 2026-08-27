"""Choosing one transcript per gene (GitHub issue #48).

The old rule kept whichever transcript was read first and silently discarded
every mutation on the others, so two MAFs holding the same mutations in a
different row order produced different rankings.
"""

import io
from collections import namedtuple

import pytest

from mutagene.io.protein_mutations_MAF import (
    read_protein_mutations_MAF,
    select_gene_transcripts,
    transcript_length,
)


def mutations(*entries):
    """Build the {sample: {(gene, transcript, change): props}} shape of the parser."""
    result = {}
    for sample, gene, transcript, change in entries:
        result.setdefault(sample, {})[(gene, transcript, change)] = {"seq5": "AAAAA"}
    return result


class TestTranscriptLength:
    def make(self, **columns):
        return namedtuple("MAF", sorted(columns))(**columns)

    def test_cdna_position_denominator_is_the_length(self):
        assert transcript_length(self.make(cdna_position="1071/2445")) == 2445

    def test_protein_position_is_converted_to_bases(self):
        assert transcript_length(self.make(protein_position="357/814")) == 814 * 3

    def test_cds_position_is_used_when_present(self):
        assert transcript_length(self.make(cds_position="12/900")) == 900

    def test_absent_column_is_not_an_error(self):
        assert transcript_length(self.make(hugo_symbol="TP53")) is None

    def test_bare_position_without_a_length_is_ignored(self):
        assert transcript_length(self.make(cdna_position="1071")) is None

    def test_unparseable_length_is_ignored(self):
        assert transcript_length(self.make(cdna_position="1071/?")) is None


class TestSelection:
    def test_longest_transcript_wins(self):
        muts = mutations(
            ("S1", "TP53", "SHORT", "R1A"),
            ("S1", "TP53", "SHORT", "R2A"),
            ("S1", "TP53", "LONG", "R3A"),
        )
        lengths = {("TP53", "SHORT"): 500, ("TP53", "LONG"): 2000}

        assert select_gene_transcripts(muts, lengths) == {"TP53": "LONG"}

    def test_length_beats_mutation_count(self):
        """Longest is what #48 asks for, even when another transcript has more rows."""
        muts = mutations(
            ("S1", "TP53", "SHORT", "R1A"),
            ("S1", "TP53", "SHORT", "R2A"),
            ("S2", "TP53", "SHORT", "R4A"),
            ("S1", "TP53", "LONG", "R3A"),
        )
        lengths = {("TP53", "SHORT"): 500, ("TP53", "LONG"): 2000}

        assert select_gene_transcripts(muts, lengths) == {"TP53": "LONG"}

    def test_falls_back_to_mutation_count_without_lengths(self):
        muts = mutations(
            ("S1", "TP53", "AAA", "R1A"),
            ("S1", "TP53", "BBB", "R2A"),
            ("S2", "TP53", "BBB", "R3A"),
        )

        assert select_gene_transcripts(muts, {}) == {"TP53": "BBB"}

    def test_a_transcript_of_unknown_length_loses_to_a_known_one(self):
        muts = mutations(
            ("S1", "TP53", "UNKNOWN", "R1A"),
            ("S1", "TP53", "UNKNOWN", "R2A"),
            ("S1", "TP53", "KNOWN", "R3A"),
        )

        assert select_gene_transcripts(muts, {("TP53", "KNOWN"): 900}) == {"TP53": "KNOWN"}

    def test_genes_are_chosen_independently(self):
        muts = mutations(
            ("S1", "TP53", "T1", "R1A"),
            ("S1", "KRAS", "K1", "G12D"),
        )

        assert select_gene_transcripts(muts, {}) == {"TP53": "T1", "KRAS": "K1"}


class TestOrderIndependence:
    """The defect itself: results depended on which row came first."""

    HEADER = (
        "Hugo_Symbol\tChromosome\tStart_Position\tVariant_Classification\t"
        "Transcript_ID\tHGVSc\tHGVSp_Short\tref_context\tTumor_Sample_Barcode\tcDNA_position"
    )

    def row(self, transcript, change="c.C4T", protein="p.R2C", length="100"):
        # ref_context is 10 + 1 + 10 bases; CGA is the arginine codon at the centre
        context = "AAAAAAAAAA" + "CGA" + "AAAAAAAA"
        return (
            f"TP53\t17\t7578406\tMissense_Mutation\t{transcript}\t{change}\t"
            f"{protein}\t{context}\tSAMPLE1\t4/{length}"
        )

    def read(self, *rows):
        text = "\n".join([self.HEADER, *rows]) + "\n"
        return read_protein_mutations_MAF(io.StringIO(text), "MAF")

    def test_same_transcript_chosen_whichever_row_is_first(self):
        long_first, _ = self.read(
            self.row("LONG", length="2000"),
            self.row("SHORT", change="c.C7T", protein="p.R3C", length="300"),
        )
        short_first, _ = self.read(
            self.row("SHORT", change="c.C7T", protein="p.R3C", length="300"),
            self.row("LONG", length="2000"),
        )

        assert {t["transcript"] for t in long_first.values()} == {"LONG"}
        assert long_first.keys() == short_first.keys()
        assert {t["transcript"] for t in short_first.values()} == {"LONG"}

    def test_the_chosen_transcript_is_reported(self):
        flat, _ = self.read(self.row("ENST00000269305", length="2000"))

        assert flat[("TP53", "R2C")]["transcript"] == "ENST00000269305"

    def test_dropped_mutations_are_counted(self):
        _flat, stats = self.read(
            self.row("LONG", length="2000"),
            self.row("SHORT", change="c.C7T", protein="p.R3C", length="300"),
        )

        assert stats["dropped_other_transcripts"] == 1


@pytest.mark.parametrize("lengths", [{}, {("TP53", "AAA"): 100, ("TP53", "BBB"): 100}])
def test_ties_are_broken_deterministically(lengths):
    muts = mutations(("S1", "TP53", "AAA", "R1A"), ("S1", "TP53", "BBB", "R2A"))

    first = select_gene_transcripts(muts, lengths)
    second = select_gene_transcripts(muts, lengths)

    assert first == second


class TestAnnotationColumnFallback:
    """A blank column must not shadow a populated one.

    MAFs commonly carry several annotation columns with only one filled in.
    Assigning each in turn let a later blank overwrite an earlier value, and the
    row was then skipped as unannotated.
    """

    HEADER = (
        "Hugo_Symbol\tChromosome\tStart_Position\tVariant_Classification\t"
        "Transcript_ID\tcDNA_Change\tHGVSc\tProtein_Change\tHGVSp_Short\t"
        "ref_context\tTumor_Sample_Barcode"
    )

    def row(self, cdna, hgvsc, protein, hgvsp_short):
        context = "AAAAAAAAAA" + "CGA" + "AAAAAAAA"
        return (
            f"TP53\t17\t7578406\tMissense_Mutation\tENST1\t{cdna}\t{hgvsc}\t"
            f"{protein}\t{hgvsp_short}\t{context}\tSAMPLE1"
        )

    def read(self, row):
        text = "\n".join([self.HEADER, row]) + "\n"
        return read_protein_mutations_MAF(io.StringIO(text), "MAF")

    def test_an_empty_later_column_does_not_shadow_an_earlier_one(self):
        flat, _stats = self.read(self.row("c.C4T", "", "p.R2C", ""))

        assert ("TP53", "R2C") in flat, "a blank HGVSc discarded a populated cDNA_Change"

    def test_a_populated_later_column_is_still_usable(self):
        flat, _stats = self.read(self.row("", "c.C4T", "", "p.R2C"))

        assert ("TP53", "R2C") in flat

    def test_whitespace_counts_as_blank(self):
        flat, _stats = self.read(self.row("c.C4T", "   ", "p.R2C", "  "))

        assert ("TP53", "R2C") in flat


def test_first_populated_prefers_the_earliest_non_blank():
    from collections import namedtuple

    from mutagene.io.protein_mutations_MAF import first_populated

    Row = namedtuple("Row", "a b c")
    assert first_populated(Row("", "second", "third"), "a", "b", "c") == "second"
    assert first_populated(Row("first", "", ""), "a", "b", "c") == "first"
    assert first_populated(Row("", "", ""), "a", "b", "c") is None
    assert first_populated(Row("", "", ""), "missing") is None
