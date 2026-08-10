"""Tests for reading cohort bundles."""

import tarfile

import pytest

from mutagene.io.cohorts import (
    list_cohorts_in_tar,
    read_cohort_mutations_from_tar,
    read_cohort_size_from_profile_str,
)


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
        # Distinct stems per file type, so that listing the wrong file type
        # changes the result instead of being hidden by deduplication.
        tar_fname = bundle(
            [
                "cohorts/LUAD.aa_mutations.txt",
                "cohorts/BRCA.profile",
                "cohorts/SKCM.dna_mutations.txt",
                "cohorts/notes.txt",
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


class TestListingMatchesLoading:
    """Every cohort that gets listed has to be loadable under that same name."""

    @pytest.mark.parametrize("prefix", ["cohorts/", "cohorts/nested/", ""])
    def test_listed_cohorts_can_be_loaded(self, tmp_path, prefix):
        tar_fname = tmp_path / "cohorts.tar.gz"
        aa = tmp_path / "aa.txt"
        aa.write_text("TP53 R175H 12\n")
        profile = tmp_path / "p.profile"
        profile.write_text("#NSAMPLES 7\n")
        with tarfile.open(tar_fname, "w:gz") as tar:
            tar.add(aa, arcname=f"{prefix}LUAD.aa_mutations.txt")
            tar.add(profile, arcname=f"{prefix}LUAD.profile")

        listed = [line.strip() for line in list_cohorts_in_tar(tar_fname).splitlines()]
        assert listed == ["LUAD"]

        _, cohort_size, aa_mutations, _ = read_cohort_mutations_from_tar(tar_fname, listed[0])
        assert cohort_size == 7
        assert aa_mutations["TP53"]["R175H"] == 12

    def test_lookup_is_case_insensitive(self, bundle):
        tar_fname = bundle(["cohorts/LUAD.aa_mutations.txt"])
        _, _, aa_mutations, _ = read_cohort_mutations_from_tar(tar_fname, "luad")
        assert aa_mutations["TP53"]["R175H"] == 12

    @pytest.mark.parametrize("requested", ["AD", "LUA", "LUADX", "LUAD_extra"])
    def test_only_the_exact_stem_matches(self, bundle, requested):
        # Neither a substring nor a prefix of a cohort name may load it: the
        # stem has to match in full.
        tar_fname = bundle(["cohorts/LUAD.aa_mutations.txt"])
        _, _, aa_mutations, _ = read_cohort_mutations_from_tar(tar_fname, requested)
        assert aa_mutations == {}

    def test_neighbouring_cohort_is_not_picked_up(self, bundle):
        # 'LUAD' and 'LUAD_TCGA' coexist; asking for one must not load the other.
        tar_fname = bundle(
            [
                "cohorts/LUAD.aa_mutations.txt",
                "cohorts/LUAD_TCGA.aa_mutations.txt",
            ]
        )
        assert list_cohorts_in_tar(tar_fname) == "\tLUAD\n\tLUAD_TCGA"
        _, _, aa_mutations, _ = read_cohort_mutations_from_tar(tar_fname, "LUAD")
        assert aa_mutations["TP53"]["R175H"] == 12


class TestCohortSize:
    def test_reads_numeric_cohort_size(self):
        assert read_cohort_size_from_profile_str("#NSAMPLES 42\nA[C>A]A 1\n") == 42

    @pytest.mark.parametrize("value", ["None", "n/a", "", "12.5"])
    def test_non_numeric_cohort_size_is_unknown(self, value):
        # Published cohorts carry '#NSAMPLES None'; an unknown size must not
        # raise and take down the whole run.
        assert read_cohort_size_from_profile_str(f"#NSAMPLES {value}\nA[C>A]A 1\n") == 0

    def test_missing_cohort_size_is_zero(self):
        assert read_cohort_size_from_profile_str("A[C>A]A 1\n") == 0

    def test_tolerates_other_comment_lines(self):
        profile_str = "# generated by mutagene on some date\n#NSAMPLES 7\n"
        assert read_cohort_size_from_profile_str(profile_str) == 7
