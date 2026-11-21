"""Convert Cancer Hotspots v2 alleles to a MutaGene reference TSV."""

import argparse
import re

import pandas as pd

from mutagene.mutability.enrichment import is_missense


def convert(source, output):
    table = pd.read_excel(source, sheet_name="SNV-hotspots", dtype=str)
    variants = set()
    for row in table.itertuples():
        position = str(row.Amino_Acid_Position)
        if not re.fullmatch(r"[1-9][0-9]*", position):
            continue
        ref = row.Reference_Amino_Acid.split(":")[0]
        for alt_counts in row.Variant_Amino_Acid.split("|"):
            change = ref + position + alt_counts.split(":")[0]
            if is_missense(change):
                variants.add((row.Hugo_Symbol, change))
    pd.DataFrame(sorted(variants), columns=["gene", "mutation"]).to_csv(
        output, sep="\t", index=False
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source")
    parser.add_argument("output")
    args = parser.parse_args()
    convert(args.source, args.output)
