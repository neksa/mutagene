"""Compare missense driver candidates with known genes and exact protein changes."""

import json
import re
from importlib.resources import files

import pandas as pd
from scipy.stats import fisher_exact

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
MISSENSE = re.compile(rf"([{AMINO_ACIDS}])[1-9][0-9]*([{AMINO_ACIDS}])")


def is_missense(change):
    match = MISSENSE.fullmatch(str(change))
    return match is not None and match[1] != match[2]


def read_known_variants(path):
    variants = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if not {"gene", "mutation"}.issubset(variants):
        raise ValueError("Known variants require gene and mutation columns")
    if (
        variants.empty
        or not variants.mutation.map(is_missense).all()
        or (variants.gene.str.strip() == "").any()
    ):
        raise ValueError("Supply gene symbols and missense protein changes such as PIK3CA E545K")
    return set(zip(variants.gene, variants.mutation))


def compare(unit, driver, other, known):
    hits = len(driver & known)
    other_hits = len(other & known)
    fold = p = None
    if driver and other:
        if other_hits:
            fold = (hits / len(driver)) / (other_hits / len(other))
        p = fisher_exact(
            [[hits, len(driver) - hits], [other_hits, len(other) - other_hits]],
            alternative="greater",
        )[1]
    return dict(
        unit=unit,
        driver_hits=hits,
        driver_total=len(driver),
        other_hits=other_hits,
        other_total=len(other),
        fold_enrichment=fold,
        fisher_p=p,
    )


def driver_enrichment(ranked, known_variants=None):
    """Count each eligible gene or gene–protein change once."""
    genes_path = files("mutagene").joinpath("data", "known_drivers.json")
    known_genes = set(json.loads(genes_path.read_text())["genes"])
    scorable = ranked.loc[
        ranked.label.isin(["Driver", "Potential driver", "Passenger"])
        & ranked.mutation.map(is_missense)
    ]
    driver = scorable.loc[scorable.label == "Driver"]
    genes = set(driver.gene)
    rows = [compare("gene", genes, set(scorable.gene) - genes, known_genes)]
    if known_variants is not None:
        alleles = set(zip(driver.gene, driver.mutation))
        other = set(zip(scorable.gene, scorable.mutation)) - alleles
        rows.append(compare("variant", alleles, other, known_variants))
    return pd.DataFrame(rows)
