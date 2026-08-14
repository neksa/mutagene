"""Regression tests for MAF header parsing (GitHub issue #91)."""

import io

import pytest

from mutagene.io import mutations_profile
from mutagene.io.context_stats import new_context_stats

MAF_COLUMNS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Strand",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
]

MAF_ROWS = [
    ["TP53", "17", "7578406", "7578406", "+", "SNP", "C", "C", "T", "SAMPLE1"],
    ["KRAS", "12", "25398284", "25398284", "+", "SNP", "C", "C", "T", "SAMPLE1"],
]


def build_maf(n_empty_trailing_columns=0, newline="\n"):
    """Build MAF text with an optional run of empty trailing columns."""
    header = MAF_COLUMNS + [""] * n_empty_trailing_columns
    lines = ["\t".join(header)]
    lines += ["\t".join(row + [""] * n_empty_trailing_columns) for row in MAF_ROWS]
    return newline.join(lines) + newline


@pytest.fixture
def stub_context(monkeypatch):
    """Return a fixed 5'/3' context so tests do not need a genome file."""

    def fake_get_context_batch(mutations, assembly, method="twobit"):
        contexts = {(chrom, pos): ("A", "G") for chrom, pos, _x, _y in mutations}
        return contexts, new_context_stats()

    monkeypatch.setattr(mutations_profile, "get_context_batch", fake_get_context_batch)


def test_maf_header_with_trailing_empty_columns(stub_context):
    """Header ending in empty tab-delimited columns must not truncate the header."""
    muts = io.StringIO(build_maf(n_empty_trailing_columns=7))
    mutations, stats = mutations_profile.read_auto_profile(muts, fmt=None, asm=None)

    assert stats["format"] == "MAF"
    assert stats["loaded"] == 2
    assert mutations["AGCT"] == 2.0


def test_normal_maf_without_trailing_columns(stub_context):
    muts = io.StringIO(build_maf())
    mutations, stats = mutations_profile.read_auto_profile(muts, fmt=None, asm=None)

    assert stats["format"] == "MAF"
    assert stats["loaded"] == 2
    assert mutations["AGCT"] == 2.0


def test_maf_with_crlf_line_endings(stub_context):
    muts = io.StringIO(build_maf(n_empty_trailing_columns=3, newline="\r\n"))
    mutations, stats = mutations_profile.read_auto_profile(muts, fmt=None, asm=None)

    assert stats["format"] == "MAF"
    assert stats["loaded"] == 2
    assert mutations["AGCT"] == 2.0


def test_maf_with_version_comment(stub_context):
    """The '#version 2.x' sniffing shortcut still finds the header below it."""
    text = "#version 2.4\n" + build_maf(n_empty_trailing_columns=4)
    muts = io.StringIO(text)
    mutations, stats = mutations_profile.read_auto_profile(muts, fmt=None, asm=None)

    assert stats["format"] == "MAF"
    assert stats["loaded"] == 2
    assert mutations["AGCT"] == 2.0


def test_crlf_rows_are_parsed_without_line_terminators(monkeypatch):
    """Data rows must not carry \\r or \\n into their fields."""
    parsed = []

    def capture_context_batch(mutations, assembly, method="twobit"):
        parsed.extend(mutations)
        contexts = {(chrom, pos): ("A", "G") for chrom, pos, _x, _y in mutations}
        return contexts, new_context_stats()

    monkeypatch.setattr(mutations_profile, "get_context_batch", capture_context_batch)

    text = "\t".join(MAF_COLUMNS) + "\r\n" + "\t".join(MAF_ROWS[0]) + "\r\n"
    mutations_profile.read_auto_profile(io.StringIO(text), fmt=None, asm=None)

    assert parsed == [("17", 7578406, "C", "T")]


def test_vcf_sniffing_still_works(stub_context):
    text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr17\t7578406\t.\tC\tT\t.\tPASS\t.\n"
        "chr12\t25398284\t.\tC\tT\t.\tPASS\t.\n"
    )
    muts = io.StringIO(text)
    mutations, stats = mutations_profile.read_auto_profile(muts, fmt=None, asm=None)

    assert stats["format"] == "VCF"
    assert stats["loaded"] == 2
    assert mutations["AGCT"] == 2.0


def test_explicit_maf_format_with_trailing_empty_columns(stub_context):
    """fmt='MAF' skips sniffing but must handle the same header shape."""
    muts = io.StringIO(build_maf(n_empty_trailing_columns=7))
    mutations, stats = mutations_profile.read_auto_profile(muts, fmt="MAF", asm=None)

    assert stats["format"] == "MAF"
    assert stats["loaded"] == 2
    assert mutations["AGCT"] == 2.0
