"""Keep annotation columns separate from ranking statistics."""

import io

import pandas as pd
import pytest

from mutagene.mutability.mutability import rank


def write_ranking(annotations=None):
    mutations = {
        (gene, "R248W"): {"seq5": {"ACGGA": 1}, "transcript": "T1"} for gene in ("FIRST", "SECOND")
    }
    stream = io.StringIO()
    rank(mutations, stream, [10.0] * 96, None, 1, 8.03e-05, 0.0034, annotations=annotations)
    stream.seek(0)
    return pd.read_csv(stream, sep="\t")


def test_later_annotation_columns_survive_without_changing_scores():
    baseline = write_ranking()
    annotated = write_ranking(
        {
            ("FIRST", "R248W"): {"alphamissense_status": "missing"},
            ("SECOND", "R248W"): {
                "alphamissense_status": "matched",
                "alphamissense_transcript": "T1",
                "alphamissense_uniprot": "P04637",
                "alphamissense_score": 0.9,
            },
        }
    )
    pd.testing.assert_frame_equal(baseline, annotated[baseline.columns])
    second = annotated.set_index("gene").loc["SECOND"]
    assert second["alphamissense_transcript"] == "T1"
    assert second["alphamissense_uniprot"] == "P04637"
    assert second["alphamissense_score"] == pytest.approx(0.9)


@pytest.mark.parametrize("column", ["bscore", "qvalue", "label", "gene"])
def test_annotations_cannot_replace_ranking_columns(column):
    with pytest.raises(ValueError, match="must not replace ranking"):
        write_ranking({("FIRST", "R248W"): {column: 0}})
