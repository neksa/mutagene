from MAFmutationmotifs import get_enrichment


def test_MAF():
    observed = get_enrichment("short_mut.txt", "TCW", 1, "C", "T", 19, assembly=37)

    expected = (0.893104332222507, 0, 1.0, 0.9328516465827816, 6, 47, 249, 1742)

    assert observed == expected
