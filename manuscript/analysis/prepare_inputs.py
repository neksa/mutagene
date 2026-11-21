"""Prepare small manuscript inputs from the BGI and SignatureAnalyzer source files."""

import argparse
import hashlib
from pathlib import Path

import pandas as pd

from mutagene.io.protein_mutations_MAF import read_protein_mutations_MAF

HERE = Path(__file__).resolve().parents[1]
RANKED = HERE.parent / "analysis/release-comparison/pcawg-penta/bgi-ranked.tsv.gz"


def prepare(bgi, genome, counts, clinical, references):
    if (
        hashlib.sha256(bgi.read_bytes()).hexdigest()
        != "c09b911056ba5d5f750df589a856706a1eb782290b78bdf32b68521fcfffe7d7"
    ):
        raise ValueError("Expected the published BGI mutation table")
    out = HERE / "inputs"
    with bgi.open() as handle:
        proteins, _ = read_protein_mutations_MAF(handle, str(genome))
    ranked = pd.read_csv(RANKED, sep="\t", comment="#")
    rows = []
    for r in ranked.itertuples():
        contexts = proteins[r.gene, r.mutation]["seq5"]
        assert sum(contexts.values()) == r.observed
        seq = min(contexts, key=lambda k: (-contexts[k], k))
        rows.append((r.gene, r.mutation, seq))
    pd.DataFrame(rows, columns=["gene", "mutation", "seq5"]).to_csv(
        out / "bgi-codons.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0}
    )
    maf = pd.read_csv(bgi, sep="\t", comment="#", low_memory=False)
    maf["mutation"] = maf.HGVSp_Short.str.removeprefix("p.")
    keep = maf.rename(
        columns={
            "Hugo_Symbol": "gene",
            "Transcript_ID": "transcript",
            "Tumor_Sample_Barcode": "sample",
        }
    ).merge(ranked[["gene", "mutation", "transcript"]], on=["gene", "mutation", "transcript"])
    carriers = keep[["gene", "mutation", "sample"]].drop_duplicates()
    check = ranked.merge(
        carriers.groupby(["gene", "mutation"]).size().rename("carriers"),
        on=["gene", "mutation"],
        how="left",
    )
    assert (check.observed == check.carriers).all()
    carriers.to_csv(
        out / "bgi-carriers.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    tumors = pd.read_csv(clinical, sep="\t", comment="#")[["SAMPLE_ID", "TNM"]]
    profiles = pd.read_csv(counts, sep="\t", index_col=0)
    assert len(tumors) == 99 and set(tumors.SAMPLE_ID) == set(profiles.columns)
    tumors.to_csv(out / "bgi-stage.tsv", sep="\t", index=False)
    profiles.to_csv(out / "bgi-sbs96.tsv", sep="\t")
    for context, cols in [("SBS1536", ["SBS39_P", "SBS17b_P", "SBS7a_S"]), ("ID83", ["ID1"])]:
        ref = pd.read_csv(references / f"{context}-reference.tsv", sep="\t", index_col=0)
        ref[cols].to_csv(out / f"{context}-display.tsv", sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("bgi", "genome", "counts", "clinical", "references"):
        parser.add_argument(name, type=Path)
    prepare(**vars(parser.parse_args()))
