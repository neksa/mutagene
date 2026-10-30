"""Variant annotation must not silently substitute alleles or transcripts."""

import pytest

from mutagene.io.alphamissense import AlphaMissenseLookup, annotate_protein_mutations

pysam = pytest.importorskip("pysam")


@pytest.fixture
def indexed_scores(tmp_path):
    table = tmp_path / "scores.tsv"
    table.write_text(
        "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n"
        "chr1\t69094\tG\tT\thg38\tQ8NH21\tENST00000335137.4\tV2L\t0.2937\tlikely_benign\n"
        "chr1\t69094\tG\tA\thg38\tQ8NH21\tENST00000335137.4\tV2M\t0.3296\tlikely_benign\n"
    )
    compressed = str(table) + ".gz"
    pysam.tabix_compress(str(table), compressed)
    pysam.tabix_index(compressed, seq_col=0, start_col=1, end_col=1, zerobased=False)
    return compressed


def test_exact_variant_returns_score(indexed_scores):
    with AlphaMissenseLookup(indexed_scores, "hg38") as lookup:
        result = lookup.lookup("1", 69094, "G", "T", "ENST00000335137.4", "V2L")
    assert result["alphamissense_score"] == pytest.approx(0.2937)
    assert result["alphamissense_status"] == "matched"


@pytest.mark.parametrize(
    "transcript,change,status",
    [
        ("ENST00000335137.5", "V2L", "transcript_mismatch"),
        ("ENST00000335137.4", "V2M", "protein_mismatch"),
        ("ENST00000335137", "V2L", "matched_unversioned_transcript"),
    ],
)
def test_mapping_status_is_explicit(indexed_scores, transcript, change, status):
    with AlphaMissenseLookup(indexed_scores, "hg38") as lookup:
        result = lookup.lookup("1", 69094, "G", "T", transcript, change)
    assert result["alphamissense_status"] == status


def test_wrong_assembly_fails(indexed_scores):
    with AlphaMissenseLookup(indexed_scores, "hg19") as lookup:
        with pytest.raises(ValueError, match="different genome"):
            lookup.lookup("1", 69094, "G", "T", "ENST00000335137.4", "V2L")


def test_missing_is_not_a_benign_score(indexed_scores):
    with AlphaMissenseLookup(indexed_scores, "hg38") as lookup:
        result = lookup.lookup("1", 69094, "G", "C", "ENST00000335137.4", "V2L")
    assert result["alphamissense_status"] == "missing"
    assert result["alphamissense_score"] is None


def test_annotation_retains_separate_evidence(indexed_scores):
    mutations = {
        ("GENE", "V2L"): {
            "transcript": "ENST00000335137.4",
            "genomic_variants": {("1", 69094, "G", "T")},
        }
    }
    with AlphaMissenseLookup(indexed_scores, "hg38") as lookup:
        annotations = annotate_protein_mutations(mutations, lookup)
    assert annotations[("GENE", "V2L")]["alphamissense_score"] == pytest.approx(0.2937)
    assert "bscore" not in annotations[("GENE", "V2L")]
    assert "label" not in annotations[("GENE", "V2L")]


@pytest.mark.parametrize("resume,corrupt", [(False, False), (True, False), (True, True)])
def test_download_verifies_bytes_before_publishing(tmp_path, monkeypatch, resume, corrupt):
    import hashlib
    import io
    from pathlib import Path

    from mutagene.io import alphamissense as module

    payload = b"verified tiny table"
    directory = tmp_path / "hg38"
    directory.mkdir()
    destination = directory / "AlphaMissense_hg38.tsv.gz"
    partial = Path(str(destination) + ".part")
    offset = 5 if resume else 0
    if resume:
        partial.write_bytes(payload[:offset])
    monkeypatch.setitem(
        module.SOURCES,
        "hg38",
        {
            "size": len(payload),
            "md5": hashlib.md5(payload).hexdigest(),
        },
    )

    class Response:
        status_code = 206
        headers = {"Content-Range": f"bytes {offset}-{len(payload)-1}/{len(payload)}"}
        raw = io.BytesIO(b"x" * (len(payload) - offset) if corrupt else payload[offset:])

        def raise_for_status(self):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

    def get(url, headers, **kwargs):
        assert headers["Range"] == f"bytes={offset}-"
        return Response()

    class Indexer:
        @staticmethod
        def tabix_index(path, **kwargs):
            Path(path + ".tbi").write_bytes(b"index")

    monkeypatch.setattr(module.requests, "get", get)
    monkeypatch.setattr(module, "require_pysam", lambda: Indexer())
    if corrupt:
        with pytest.raises(ValueError, match="checksum"):
            module.fetch_alphamissense("hg38", tmp_path)
        assert not destination.exists()
        assert not Path(str(destination) + ".tbi").exists()
    else:
        assert module.fetch_alphamissense("hg38", tmp_path) == destination
        assert destination.read_bytes() == payload
        assert Path(str(destination) + ".provenance.json").exists()


def test_invalid_manifest_does_not_open_the_tabix_file(indexed_scores, monkeypatch):
    import json
    from pathlib import Path

    Path(indexed_scores + ".provenance.json").write_text(json.dumps({"genome": "hg19"}))
    opened = []
    monkeypatch.setattr(pysam, "TabixFile", lambda path: opened.append(path))
    with pytest.raises(ValueError, match="requested assembly"):
        AlphaMissenseLookup(indexed_scores, "hg38")
    assert opened == []


def test_malformed_record_returns_a_validation_error(tmp_path):
    table = tmp_path / "malformed.tsv"
    table.write_text("chr1\t69094\tG\tT\thg38\n")
    compressed = str(table) + ".gz"
    pysam.tabix_compress(str(table), compressed)
    pysam.tabix_index(compressed, seq_col=0, start_col=1, end_col=1, zerobased=False)
    with AlphaMissenseLookup(compressed, "hg38") as lookup:
        with pytest.raises(ValueError, match="number of columns"):
            lookup.lookup("1", 69094, "G", "T", "ENST00000335137.4", "V2L")
