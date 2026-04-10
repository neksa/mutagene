"""Tests for webapp analysis helpers."""

import gzip
import tarfile
import tempfile
from pathlib import Path

import pytest

from mutagene.webapp.analysis import extract_input_file, open_input_file


class TestOpenInputFile:
    def test_plain_text(self, tmp_path):
        f = tmp_path / "input.maf"
        f.write_text("Hugo_Symbol\tChromosome\n")
        with open_input_file(f, "rt") as fh:
            assert "Hugo_Symbol" in fh.read()

    def test_gzipped(self, tmp_path):
        f = tmp_path / "input.maf.gz"
        with gzip.open(f, "wt") as gz:
            gz.write("Hugo_Symbol\tChromosome\n")
        with open_input_file(f, "rt") as fh:
            assert "Hugo_Symbol" in fh.read()


class TestExtractInputFile:
    def test_plain_file_passthrough(self, tmp_path):
        f = tmp_path / "input.maf"
        f.write_text("data")
        result = extract_input_file(f, tmp_path / "output")
        assert result == f

    def test_tar_gz_extraction(self, tmp_path):
        # Create a tar.gz with a mutation file inside
        maf_content = "Hugo_Symbol\tChromosome\nTP53\t17\n"
        maf_path = tmp_path / "mutations.maf"
        maf_path.write_text(maf_content)

        tar_path = tmp_path / "dataset.tar.gz"
        with tarfile.open(tar_path, "w:gz") as tar:
            tar.add(str(maf_path), arcname="data_mutations.maf")

        output_dir = tmp_path / "output"
        output_dir.mkdir()
        result = extract_input_file(tar_path, output_dir)
        assert result.exists()
        assert "Hugo_Symbol" in result.read_text()
        # Extracted into output dir, not /tmp
        assert str(output_dir) in str(result)

    def test_tar_gz_no_maf(self, tmp_path):
        # Create tar.gz with no mutation files
        txt = tmp_path / "readme.txt"
        txt.write_text("not a mutation file")
        tar_path = tmp_path / "archive.tar.gz"
        with tarfile.open(tar_path, "w:gz") as tar:
            tar.add(str(txt), arcname="readme.txt")

        with pytest.raises(ValueError, match="No mutation file"):
            extract_input_file(tar_path, tmp_path / "output")

    def test_tar_gz_path_traversal_stripped(self, tmp_path):
        # Tar member with path traversal gets basename-stripped, staying inside output dir
        maf_content = "Hugo_Symbol\n"
        maf_path = tmp_path / "mutations.maf"
        maf_path.write_text(maf_content)

        tar_path = tmp_path / "evil.tar.gz"
        with tarfile.open(tar_path, "w:gz") as tar:
            tar.add(str(maf_path), arcname="../../etc/mutation_passwd.maf")

        output_dir = tmp_path / "output"
        output_dir.mkdir()
        result = extract_input_file(tar_path, output_dir)
        # Should extract safely inside output dir, not to ../../etc/
        assert str(output_dir) in str(result)
        assert "etc" not in str(result)
