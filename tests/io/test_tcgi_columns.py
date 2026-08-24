"""TCGI column names (the documented CHR was never accepted).

The format documentation names CHR, POS, REF and ALT as the mandatory columns,
but the reader only ever looked for `chrom`. A file written to the documentation
was rejected with "Chromosome is not defined in TCGI file", and the column the
code actually wanted was documented nowhere.
"""

import pytest

from mutagene.io import context_window
from mutagene.io.context_stats import new_context_stats


@pytest.fixture
def stub_genome(monkeypatch):
    def get_context(mutations, twobit_file, window_size):
        contexts = {
            (chrom, pos): (("A", "G"), [(chrom, pos, x, strand)])
            for chrom, pos, strand, x, y in mutations
        }
        return contexts, new_context_stats()

    monkeypatch.setattr(context_window, "get_context_twobit_window", get_context)


def tcgi(chrom_column, extra_columns="", extra_values=""):
    header = f"{chrom_column}\tPOS\tREF\tALT{extra_columns}"
    row = f"17\t7578406\tC\tT{extra_values}"
    return [header + "\n", row + "\n"]


def read(lines):
    return context_window.read_TCGI_with_context_window(lines, "unused.2bit", 1)


@pytest.mark.parametrize("column", ["CHR", "CHROM", "chr", "chrom"])
def test_both_documented_and_actual_names_are_accepted(stub_genome, column):
    mutations, _, stats = read(tcgi(column))

    assert stats["loaded"] == 1
    assert stats["format"] == "TCGI"
    assert mutations["TCGI"]["AGCT"] == 1.0


def test_sample_column_is_still_optional(stub_genome):
    mutations, _, stats = read(tcgi("CHR"))

    assert list(mutations) == ["TCGI"], "without SAMPLE everything is one sample"


def test_sample_column_is_used_when_present(stub_genome):
    mutations, _, stats = read(tcgi("CHR", "\tSAMPLE", "\tPATIENT1"))

    assert list(mutations) == ["PATIENT1"]


def test_missing_chromosome_column_names_the_documented_one(stub_genome):
    lines = ["POS\tREF\tALT\n", "7578406\tC\tT\n"]

    with pytest.raises(ValueError, match="CHR"):
        read(lines)
