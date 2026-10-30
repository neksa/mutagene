"""VCF readers must preserve sample meaning and reject unsupported cohorts."""

import io

import pytest

from mutagene.io.context_window import read_mutations
from mutagene.io.mutations_profile import read_auto_profile
from tests.cli.local_fixtures import LOCAL_GENOME


@pytest.mark.parametrize("reader", ["profile", "context"])
def test_rejects_multiple_samples_before_context_lookup(reader):
    text = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
    with pytest.raises(ValueError, match="Multiple VCF samples"):
        if reader == "profile":
            read_auto_profile(io.StringIO(text), "VCF", LOCAL_GENOME)
        else:
            read_mutations("VCF", io.StringIO(text), LOCAL_GENOME, window_size=1)


@pytest.mark.parametrize("reader", ["profile", "context"])
def test_reference_genotype_does_not_count_as_a_mutation(reader):
    text = (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n"
        "chr1\t100\t.\tC\tT\t.\tPASS\t.\tGT\t0/0\n"
    )
    if reader == "profile":
        _, stats = read_auto_profile(io.StringIO(text), "VCF", LOCAL_GENOME)
    else:
        _, _, stats = read_mutations("VCF", io.StringIO(text), LOCAL_GENOME, window_size=1)
    assert stats["loaded"] == 0
