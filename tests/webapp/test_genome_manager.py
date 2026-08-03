"""Tests for webapp genome reference management."""

import pytest

from mutagene.webapp.genome_manager import GenomeManager


@pytest.fixture
def manager(tmp_path):
    return GenomeManager(genomes_dir=tmp_path / "genomes")


class TestGenomePathValidation:
    def test_supported_genome_resolves_inside_dir(self, manager):
        path = manager.get_genome_path("hg19")
        assert path.parent == manager.genomes_dir
        assert path.name == "hg19.2bit"

    @pytest.mark.parametrize(
        "genome",
        ["../../etc/passwd", "/etc/passwd", "hg19/../../secret", "", "HG19"],
    )
    def test_rejects_unsupported_names(self, manager, genome):
        # An unchecked name here becomes part of a filesystem path.
        with pytest.raises(ValueError):
            manager.get_genome_path(genome)

    def test_is_downloaded_is_false_for_unknown_name(self, manager):
        assert manager.is_downloaded("../../etc/passwd") is False


class TestAtomicDownload:
    def test_failed_download_leaves_no_usable_file(self, manager, monkeypatch):
        def boom(url, dst):
            # Simulate a download that dies partway through.
            with open(dst, "wb") as fh:
                fh.write(b"partial")
            raise OSError("connection reset")

        monkeypatch.setattr("mutagene.webapp.genome_manager.download_from_url", boom)

        assert manager.download_genome("hg19") is False
        assert not manager.is_downloaded("hg19")
        assert not manager.get_genome_path("hg19").exists()
        assert list(manager.genomes_dir.iterdir()) == []

    def test_successful_download_is_promoted(self, manager, monkeypatch):
        def ok(url, dst):
            with open(dst, "wb") as fh:
                fh.write(b"2bit-data")

        monkeypatch.setattr("mutagene.webapp.genome_manager.download_from_url", ok)

        assert manager.download_genome("hg19") is True
        assert manager.is_downloaded("hg19")
        assert manager.get_genome_path("hg19").read_bytes() == b"2bit-data"

    def test_download_writes_to_temp_path_first(self, manager, monkeypatch):
        seen = {}

        def capture(url, dst):
            seen["dst"] = dst
            with open(dst, "wb") as fh:
                fh.write(b"x")

        monkeypatch.setattr("mutagene.webapp.genome_manager.download_from_url", capture)
        manager.download_genome("hg38")

        assert seen["dst"].endswith(".part")
        assert seen["dst"] != str(manager.get_genome_path("hg38"))

    def test_unsupported_genome_is_not_downloaded(self, manager):
        assert manager.download_genome("nonexistent") is False


class TestGenomeInfo:
    def test_reports_all_supported_genomes(self, manager):
        info = manager.get_genome_info()
        assert set(info) == set(GenomeManager.SUPPORTED_GENOMES)
        assert all(entry["downloaded"] is False for entry in info.values())
