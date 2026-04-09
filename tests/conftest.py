"""Shared test fixtures for mutagene tests."""

import tempfile
from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def temp_dir():
    """Provide a temporary directory that is cleaned up after the test."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_maf_lines():
    """Minimal MAF-format lines for testing parsers."""
    header = (
        "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\t"
        "Start_Position\tEnd_Position\tStrand\tVariant_Classification\t"
        "Variant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\t"
        "dbSNP_RS\tdbSNP_Val_Status\tTumor_Sample_Barcode"
    )
    row = (
        "TP53\t7157\tBCGSC\tGRCh37\t17\t"
        "7578406\t7578406\t+\tMissense_Mutation\t"
        "SNP\tC\tC\tT\t"
        "novel\t\tSAMPLE1"
    )
    return [header + "\n", row + "\n"]


@pytest.fixture
def sample_profile_96():
    """96-channel mutational profile (all ones for simple testing)."""
    return [1.0] * 96
