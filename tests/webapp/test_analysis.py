"""Tests for webapp analysis helpers."""

import gzip
import tarfile

import pytest

from mutagene.webapp.analysis import (
    extract_input_file,
    open_input_file,
    profile_channel_order,
    profile_dict_to_array,
)


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
        escape_target = tmp_path.parent / "etc" / "mutation_passwd.maf"
        result = extract_input_file(tar_path, output_dir)

        # Should extract safely inside output dir, not to ../../etc/
        assert result.resolve().is_relative_to(output_dir.resolve())
        # The bytes must actually land at the returned path, and nowhere else.
        assert result.exists()
        assert result.read_text() == maf_content
        assert not escape_target.exists()

    def test_tar_gz_nested_member_is_readable(self, tmp_path):
        # A member inside a subdirectory must be returned at a path that exists;
        # previously the return value pointed at a basename that was never written.
        maf_path = tmp_path / "mutations.maf"
        maf_path.write_text("Hugo_Symbol\n")

        tar_path = tmp_path / "nested.tar.gz"
        with tarfile.open(tar_path, "w:gz") as tar:
            tar.add(str(maf_path), arcname="data/inner/mutations.maf")

        result = extract_input_file(tar_path, tmp_path / "output")
        assert result.exists()
        assert result.read_text() == "Hugo_Symbol\n"

    def test_tar_gz_oversized_member_rejected(self, tmp_path, monkeypatch):
        import mutagene.webapp.analysis as analysis_mod

        maf_path = tmp_path / "mutations.maf"
        maf_path.write_text("Hugo_Symbol\n" * 100)

        tar_path = tmp_path / "big.tar.gz"
        with tarfile.open(tar_path, "w:gz") as tar:
            tar.add(str(maf_path), arcname="mutations.maf")

        monkeypatch.setattr(analysis_mod, "MAX_EXTRACTED_BYTES", 10)
        with pytest.raises(ValueError, match="too large|exceeds"):
            extract_input_file(tar_path, tmp_path / "output")


class TestProfileChannelOrder:
    def test_matches_signature_matrix_row_order(self):
        from mutagene.io.profile import get_profile_attributes_dict, read_signatures

        order = profile_channel_order()
        assert len(order) == 96
        assert len(set(order)) == 96
        # Must equal the order used to build W, not alphabetical order.
        assert order == [
            f"{a['context'][0]}[{a['mutation'][0]}>{a['mutation'][1]}]{a['context'][1]}"
            for a in get_profile_attributes_dict()
        ]
        assert order != sorted(order), "alphabetical order would permute 88 of 96 channels"

        W, _ = read_signatures("COSMICv3", only=None)
        assert W.shape[0] == len(order)

    def test_dict_to_array_places_counts_in_signature_order(self):
        order = profile_channel_order()
        # Put a marker count in a channel whose alphabetical index differs.
        target = "A[C>A]T"
        arr = profile_dict_to_array({target: 7})
        assert arr.shape == (96,)
        assert arr[order.index(target)] == 7
        assert arr.sum() == 7

    def test_dict_to_array_fills_missing_channels(self):
        arr = profile_dict_to_array({"A[C>A]A": 3})
        assert arr.shape == (96,)
        assert arr.sum() == 3
