from broad_MAF_set import get_enrichment


def test_MAF():
    observed = get_enrichment("tcga_A0C8.maf", "CC", 1, "C", "T", 29, assembly=37)

    expected = (0.893104332222507, 0, 1.0, 0.9328516465827816, 6, 47, 249, 1742)

    assert observed == expected
