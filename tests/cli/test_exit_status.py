"""A run that read nothing must not report success.

All four analysis subcommands wrote an empty output file and exited 0 when the
input could not be read, so a shell pipeline carried on with an empty profile as
though the step had worked.
"""

import pytest

from tests.cli import cli_test_utils

GENOME = f"{cli_test_utils.TEST_DIR}/hg19.2bit"
SAMPLE = f"{cli_test_utils.TEST_DIR}/sample1.maf"


@pytest.fixture
def unreadable(tmp_path):
    path = tmp_path / "junk.maf"
    path.write_text("this is not a mutation file\njust prose\n")
    return str(path)


@pytest.fixture
def no_usable_chromosome(tmp_path):
    """Well-formed MAF whose chromosome is absent from the assembly."""
    path = tmp_path / "nochrom.maf"
    path.write_text(
        "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
        "Tumor_Seq_Allele2\tTumor_Sample_Barcode\n"
        "99\t100\t100\tC\tT\tSAMPLE1\n"
    )
    return str(path)


def run(command, args):
    with pytest.raises(SystemExit) as excinfo:
        cli_test_utils.run_with_args(command, args)
    return excinfo.value.code


class TestUnreadableInput:
    @pytest.mark.parametrize(
        "command,extra",
        [
            ("profile", []),
            ("signature", ["-s", "MGA"]),
            ("rank", []),
            ("motif", []),
        ],
    )
    def test_exit_is_non_zero(self, command, extra, unreadable, tmp_path, test_data):
        code = run(
            command,
            ["-i", unreadable, "-g", GENOME, "-o", str(tmp_path / "out.tsv")] + extra,
        )

        assert code not in (0, None), f"{command} reported success on an unreadable file"

    def test_a_file_with_no_usable_chromosome_also_fails(
        self, no_usable_chromosome, tmp_path, test_data
    ):
        code = run(
            "profile",
            ["-i", no_usable_chromosome, "-g", GENOME, "-o", str(tmp_path / "out.tsv")],
        )

        assert code not in (0, None)


class TestUsableInput:
    """The converse: a real file must still succeed."""

    def test_profile_succeeds(self, tmp_path, test_data):
        out = tmp_path / "profile.tsv"
        cli_test_utils.run_with_args("profile", ["-i", SAMPLE, "-g", GENOME, "-o", str(out)])

        assert out.stat().st_size > 0

    def test_signature_succeeds(self, tmp_path, test_data):
        out = tmp_path / "sig.tsv"
        cli_test_utils.run_with_args(
            "signature", ["-i", SAMPLE, "-g", GENOME, "-s", "MGA", "-o", str(out)]
        )

        assert out.stat().st_size > 0
