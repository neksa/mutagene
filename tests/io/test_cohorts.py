"""Tests for reading cohort bundles."""

import tarfile

import pytest

from mutagene.io.cohorts import list_cohorts_in_tar


def make_tar(path, names):
    """Builds a cohorts tarball containing the given archive member names."""
    with tarfile.open(path, "w:gz") as tar:
        for name in names:
            payload = path.parent / "payload.txt"
            payload.write_text("TP53 R175H 12\n")
            tar.add(payload, arcname=name)
    return path


@pytest.fixture
def bundle(tmp_path):
    def build(names):
        return make_tar(tmp_path / "cohorts.tar.gz", names)

    return build


class TestListCohortsInTar:
    def test_lists_cohorts_sorted(self, bundle):
        tar_fname = bundle(
            [
                "cohorts/LUAD.aa_mutations.txt",
                "cohorts/BRCA.aa_mutations.txt",
                "cohorts/gcb_lymphomas.aa_mutations.txt",
            ]
        )
        assert list_cohorts_in_tar(tar_fname) == "\tBRCA\n\tLUAD\n\tgcb_lymphomas"

    def test_preserves_case(self, bundle):
        # Cohort lookup is case insensitive, but the listing should show the
        # name as it is actually stored in the bundle.
        tar_fname = bundle(["cohorts/LUAD.aa_mutations.txt"])
        assert "LUAD" in list_cohorts_in_tar(tar_fname)

    def test_keeps_dots_in_cohort_name(self, bundle):
        tar_fname = bundle(["cohorts/pan.cancer.aa_mutations.txt"])
        assert list_cohorts_in_tar(tar_fname) == "\tpan.cancer"

    def test_ignores_other_cohort_files(self, bundle):
        tar_fname = bundle(
            [
                "cohorts/LUAD.aa_mutations.txt",
                "cohorts/LUAD.profile",
                "cohorts/LUAD.dna_mutations.txt",
            ]
        )
        assert list_cohorts_in_tar(tar_fname) == "\tLUAD"

    def test_skips_appledouble_sidecars(self, bundle):
        # Bundles repacked on macOS carry ._ sidecars that match the suffix.
        tar_fname = bundle(
            [
                "cohorts/LUAD.aa_mutations.txt",
                "cohorts/._LUAD.aa_mutations.txt",
            ]
        )
        assert list_cohorts_in_tar(tar_fname) == "\tLUAD"

    def test_handles_flat_bundle_without_directory(self, bundle):
        tar_fname = bundle(["LUAD.aa_mutations.txt"])
        assert list_cohorts_in_tar(tar_fname) == "\tLUAD"

    def test_empty_bundle_returns_empty_string(self, bundle):
        tar_fname = bundle(["cohorts/README.txt"])
        assert list_cohorts_in_tar(tar_fname) == ""
