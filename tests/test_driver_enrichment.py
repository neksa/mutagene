import pandas as pd
import pytest

from mutagene.mutability.enrichment import driver_enrichment, read_known_variants


def test_enrichment_uses_distinct_scorable_missense_changes():
    ranked = pd.DataFrame(
        [
            ("PIK3CA", "E545K", "Driver"),
            ("PIK3CA", "E545K", "Driver"),
            ("PIK3CA", "E545A", "Potential driver"),
            ("OTHER", "Y2A", "Passenger"),
            ("KRAS", "G12=", "Driver"),
            ("TP53", "R2R", "Driver"),
            ("RB1", "A3C", "Undefined"),
        ],
        columns=["gene", "mutation", "label"],
    )
    result = driver_enrichment(ranked, {("PIK3CA", "E545K")}).set_index("unit")
    assert list(result.loc["gene"].iloc[:4]) == [1, 1, 0, 1]
    assert list(result.loc["variant"].iloc[:4]) == [1, 1, 0, 2]
    assert pd.isna(result.loc["variant", "fold_enrichment"])
    assert result.loc["variant", "fisher_p"] == pytest.approx(1 / 3)
    driver_only = driver_enrichment(ranked.iloc[:1]).iloc[0]
    assert driver_only.other_total == 0
    assert pd.isna(driver_only.fisher_p)


def test_known_variants_require_missense_changes(tmp_path):
    reference = tmp_path / "known.tsv"
    reference.write_text("gene\tmutation\nKRAS\tG12D\n")
    assert read_known_variants(reference) == {("KRAS", "G12D")}
    reference.write_text("gene\tmutation\nKRAS\tG12G\n")
    with pytest.raises(ValueError, match="missense"):
        read_known_variants(reference)
