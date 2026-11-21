"""Check channel identity and exposure recovery for extended signature fitting."""

import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mutagene.io.matrix import channel_labels, read_count_matrix, read_reference_matrix
from mutagene.profiles.profile import read_sample_profiles
from tests.cli.cli_test_utils import run_with_args

DEMO = Path(__file__).resolve().parents[2] / "mutagene/data/demo"


@pytest.mark.parametrize(
    "context,first,second",
    [
        ("SBS1536", "AA[C>T]AA", "TA[C>T]AT"),
        ("ID83", "2:Del:R:1", "2:Del:M:1"),
    ],
)
def test_mixture_recovers_distinct_channels(tmp_path, context, first, second):
    labels = channel_labels(context)
    assert len(set(labels)) == {"SBS1536": 1536, "ID83": 83}[context]
    catalog = pd.DataFrame(0.0, index=labels, columns=["one", "two"])
    catalog.loc[first, "one"] = 1
    catalog.loc[second, "two"] = 1
    counts = pd.DataFrame({"A": 70 * catalog.one + 30 * catalog.two, "B": 100 * catalog.two})
    catalog.iloc[::-1].to_csv(tmp_path / "reference.tsv", sep="\t", index_label="Type")
    counts.to_csv(tmp_path / "counts.tsv", sep="\t", index_label="MutationType")
    run_with_args(
        "signature",
        [
            "-f",
            "MATRIX",
            "--context",
            context,
            "--reference-matrix",
            str(tmp_path / "reference.tsv"),
            "-i",
            str(tmp_path / "counts.tsv"),
            "-o",
            str(tmp_path / "fits.tsv"),
            "--diagnostics-out",
            str(tmp_path / "fit.json"),
        ],
    )
    fit = pd.read_csv(tmp_path / "fits.tsv", sep="\t").set_index(["sample", "signature"])
    assert fit.loc[("A", "one"), "exposure"] == pytest.approx(0.7)
    assert fit.loc[("A", "two"), "mutations"] == 30
    assert fit.loc[("B", "two"), "exposure"] == pytest.approx(1)
    diagnostics = json.loads((tmp_path / "fit.json").read_text())
    assert diagnostics["context"] == context
    assert diagnostics["samples"]["A"]["null_components"] == []


def test_pentanucleotide_profile_collapses_to_sbs96():
    def profile(context):
        with (DEMO / "sample.maf").open() as source:
            return read_sample_profiles(
                "MAF",
                source,
                str(DEMO / "reference.2bit"),
                context,
                keep_filtered=False,
                filter_column="FILTER",
            )

    expected = profile("SBS96")
    for sample, values in profile("SBS1536").items():
        by_triplet = (
            pd.Series(values, index=[x[1:-1] for x in channel_labels("SBS1536")])
            .groupby(level=0)
            .sum()
        )
        np.testing.assert_array_equal(by_triplet.loc[channel_labels()], expected[sample])


def test_extended_reference_rejects_wrong_scheme_and_invalid_values():
    data = pd.DataFrame({"first": 1}, index=channel_labels("ID83"))
    with pytest.raises(ValueError, match="1536"):
        read_reference_matrix(io.StringIO(data.to_csv(sep="\t", index_label="Type")), "SBS1536")
    data.iloc[0] = -1
    with pytest.raises(ValueError, match="nonnegative"):
        read_reference_matrix(io.StringIO(data.to_csv(sep="\t", index_label="Type")), "ID83")
    with pytest.raises(ValueError, match="integer counts"):
        read_count_matrix(io.StringIO(data.to_csv(sep="\t", index_label="Type")), "ID83")


def test_pentanucleotide_strand_orientation(monkeypatch):
    from mutagene.profiles import profile

    def window(sequence):
        return [("1", i, base, "+") for i, base in enumerate(sequence)]

    # All three calls describe the same pyrimidine-oriented substitution.
    records = [
        ("1", 3, "+", "C", "T", window("AGCTA")),
        ("1", 9, "-", "G", "A", window("TAGCT")),
        ("1", 15, "-", "C", "T", window("TAGCT")),
    ]
    monkeypatch.setattr(profile, "read_mutations", lambda *args, **kwargs: ({}, {"S": records}, {}))
    result = read_sample_profiles(
        "MAF", None, "", "SBS1536", keep_filtered=False, filter_column="FILTER"
    )
    assert result["S"].sum() == result["S"][channel_labels("SBS1536").index("AG[C>T]TA")] == 3


def test_variant_input_matches_matrix_fit_and_seeded_bootstrap(tmp_path):
    common = ["--context", "SBS1536", "-g", str(DEMO / "reference.2bit")]
    matrix = tmp_path / "penta.tsv"
    run_with_args(
        "profile", [*common, "--per-sample", "-i", str(DEMO / "sample.maf"), "-o", str(matrix)]
    )
    counts = pd.read_csv(matrix, sep="\t", index_col=0)
    ref = tmp_path / "catalog.tsv"
    counts.to_csv(ref, sep="\t")
    options = [*common, "--reference-matrix", str(ref), "--bootstrap", "-br", "10", "--seed", "42"]
    for name, fmt, source in [("direct", "MAF", DEMO / "sample.maf"), ("matrix", "MATRIX", matrix)]:
        run_with_args(
            "signature",
            [*options, "-f", fmt, "-i", str(source), "-o", str(tmp_path / (name + ".tsv"))],
        )
    assert (tmp_path / "direct.tsv").read_text() == (tmp_path / "matrix.tsv").read_text()
