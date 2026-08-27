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


class TestGenomeMismatchWarning:
    """The warning that could never fire (GitHub issue #99).

    The webapp counted log records around `calc_profile`, but the message it
    counted is emitted by a different function, after the handler was removed.
    Analyses against the wrong assembly therefore completed silently: the
    stored evidence was a breast-cancer cohort analysed against hg19 that
    returned haloalkane and azathioprine signatures without a word of warning.
    """

    HEADER = (
        "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
        "Tumor_Seq_Allele2\tTumor_Sample_Barcode"
    )

    def maf(self, tmp_path, ref, count=60):
        rows = [self.HEADER]
        for i in range(count):
            pos = 10 + i
            rows.append(f"17\t{pos}\t{pos}\t{ref}\tT\tSAMPLE1")
        path = tmp_path / "input.maf"
        path.write_text("\n".join(rows) + "\n")
        return path

    @pytest.fixture
    def wired_up(self, tmp_path, monkeypatch):
        """An all-A chr17 plus a genome path that exists."""
        import twobitreader

        from mutagene.webapp import analysis as analysis_module

        class FakeChromosome:
            def __init__(self, sequence):
                self.sequence = sequence

            def __getitem__(self, key):
                return self.sequence[key]

        class FakeTwoBitFile:
            def __contains__(self, name):
                return name == "chr17"

            def __getitem__(self, name):
                return FakeChromosome("A" * 200)

        monkeypatch.setattr(twobitreader, "TwoBitFile", lambda fname: FakeTwoBitFile())

        genome_path = tmp_path / "hg19.2bit"
        genome_path.write_bytes(b"stub")
        monkeypatch.setattr(
            analysis_module.GenomeManager,
            "get_genome_path",
            lambda self, genome: genome_path,
        )

    def test_wrong_assembly_is_reported(self, tmp_path, wired_up):
        from mutagene.webapp.analysis import run_cohort_analysis

        results = run_cohort_analysis(self.maf(tmp_path, ref="G"), tmp_path / "out", genome="hg19")

        warning = results.get("genome_warning")
        assert warning is not None, "a wrong assembly must not complete silently"
        assert warning["mismatch_count"] == 60
        assert "hg38" in warning["message"]

    def test_matching_assembly_is_not_reported(self, tmp_path, wired_up):
        from mutagene.webapp.analysis import run_cohort_analysis

        results = run_cohort_analysis(self.maf(tmp_path, ref="A"), tmp_path / "out", genome="hg19")

        assert "genome_warning" not in results


class TestTruncatedUploads:
    """A part-uploaded archive must say so, not leak an EOFError.

    Analysis 13 in the stored database failed with "Compressed file ended before
    the end-of-stream marker was reached", which names neither the file nor
    anything the uploader could do about it.
    """

    def gzipped(self, tmp_path, keep_bytes=None):
        import gzip as gz

        source = tmp_path / "sample.maf"
        source.write_text("Chromosome\tStart_Position\n" + "17\t100\n" * 500)
        blob = gz.compress(source.read_bytes())
        path = tmp_path / "sample.maf.gz"
        path.write_bytes(blob if keep_bytes is None else blob[:keep_bytes])
        return path

    def tarred(self, tmp_path, keep_bytes=None):
        import tarfile

        source = tmp_path / "data_mutations.maf"
        source.write_text("Chromosome\tStart_Position\n" + "17\t100\n" * 500)
        full = tmp_path / "full.tar.gz"
        with tarfile.open(full, "w:gz") as tar:
            tar.add(str(source), arcname="data_mutations.maf")
        if keep_bytes is None:
            return full
        path = tmp_path / "truncated.tar.gz"
        path.write_bytes(full.read_bytes()[:keep_bytes])
        return path

    def test_truncated_gzip_is_reported_usefully(self, tmp_path):
        from mutagene.webapp.analysis import TruncatedInputError, open_input_file

        path = self.gzipped(tmp_path, keep_bytes=40)

        with pytest.raises(TruncatedInputError) as excinfo:
            with open_input_file(path, "rt") as fh:
                for _ in fh:
                    pass

        assert "sample.maf.gz" in str(excinfo.value)
        assert "uploading it again" in str(excinfo.value)

    def test_intact_gzip_still_reads(self, tmp_path):
        from mutagene.webapp.analysis import open_input_file

        with open_input_file(self.gzipped(tmp_path), "rt") as fh:
            lines = list(fh)

        assert lines[0].startswith("Chromosome")
        assert len(lines) == 501

    def test_truncated_tarball_is_reported_usefully(self, tmp_path):
        from mutagene.webapp.analysis import TruncatedInputError, extract_input_file

        path = self.tarred(tmp_path, keep_bytes=200)

        with pytest.raises(TruncatedInputError) as excinfo:
            extract_input_file(path, tmp_path / "out")

        assert "truncated.tar.gz" in str(excinfo.value)

    def test_intact_tarball_still_extracts(self, tmp_path):
        from mutagene.webapp.analysis import extract_input_file

        extracted = extract_input_file(self.tarred(tmp_path), tmp_path / "out")

        assert extracted.read_text().startswith("Chromosome")

    def test_a_truncated_error_is_a_value_error(self, tmp_path):
        """run_analysis reports ValueError to the client; EOFError it would not."""
        from mutagene.webapp.analysis import TruncatedInputError

        assert issubclass(TruncatedInputError, ValueError)


class TestInputFormatDetection:
    """A VCF upload was profiled correctly and then read again as a MAF."""

    def vcf(self, tmp_path, name="sample.vcf", rows=40):
        lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
        for i in range(rows):
            lines.append(f"chr17\t{20 + i * 4}\t.\tC\tT\t.\tPASS\t.")
        path = tmp_path / name
        path.write_text("\n".join(lines) + "\n")
        return path

    def maf(self, tmp_path, name="sample.maf"):
        path = tmp_path / name
        path.write_text(
            "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
            "Tumor_Seq_Allele2\tTumor_Sample_Barcode\n17\t100\t100\tC\tT\tSAMPLE1\n"
        )
        return path

    def test_vcf_by_extension(self, tmp_path):
        from mutagene.webapp.analysis import detect_input_format

        assert detect_input_format(self.vcf(tmp_path)) == "VCF"

    def test_maf_by_extension(self, tmp_path):
        from mutagene.webapp.analysis import detect_input_format

        assert detect_input_format(self.maf(tmp_path)) == "MAF"

    def test_gzipped_extension_is_seen_through(self, tmp_path):
        from mutagene.webapp.analysis import detect_input_format

        assert detect_input_format(tmp_path / "sample.vcf.gz") == "VCF"

    def test_txt_upload_is_decided_by_content(self, tmp_path):
        """.txt and .tsv are in the allow-list, so the name settles nothing."""
        from mutagene.webapp.analysis import detect_input_format

        assert detect_input_format(self.vcf(tmp_path, name="mutations.txt")) == "VCF"
        assert detect_input_format(self.maf(tmp_path, name="mutations.txt")) == "MAF"

    def test_a_vcf_upload_completes(self, tmp_path, monkeypatch):
        """This raised 'Chromosome is not defined in MAF file'."""
        import twobitreader

        from mutagene.webapp import analysis as analysis_module
        from mutagene.webapp.analysis import run_cohort_analysis

        class FakeChromosome:
            def __getitem__(self, key):
                # every position reads C, matching the REF in the fixture
                length = key.stop - key.start if isinstance(key, slice) else 1
                return "C" * length

        class FakeTwoBitFile:
            def __contains__(self, name):
                return name == "chr17"

            def __getitem__(self, name):
                return FakeChromosome()

        genome = tmp_path / "hg19.2bit"
        genome.write_bytes(b"stub")
        monkeypatch.setattr(twobitreader, "TwoBitFile", lambda fname: FakeTwoBitFile())
        monkeypatch.setattr(
            analysis_module.GenomeManager, "get_genome_path", lambda self, g: genome
        )

        results = run_cohort_analysis(self.vcf(tmp_path), tmp_path / "out", genome="hg19")

        assert results["mutations"] > 0, "the VCF produced no mutations"
        assert results["samples"] >= 1
