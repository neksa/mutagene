"""Recording and replaying run parameters (GitHub issue #45)."""

import argparse
import json
import sys

import pytest

from mutagene.cli.run_params import (
    add_run_parameter_arguments,
    collect_run_parameters,
    convert_recorded_lists,
    load_run_parameters,
    peek_params_in,
    write_run_parameters,
)


def build_parser():
    """A parser shaped like the real subcommands: FileType, nargs='*', ints."""
    parser = argparse.ArgumentParser()
    add_run_parameter_arguments(parser)
    parser.add_argument("--infile", "-i", nargs="*", type=argparse.FileType("r"))
    parser.add_argument(
        "--outfile", "-o", nargs="?", type=argparse.FileType("w"), default=sys.stdout
    )
    parser.add_argument("--genome", "-g", type=str, default="hg19")
    parser.add_argument("--nsamples", "-n", type=int)
    parser.add_argument("--verbose", action="count", default=0)
    return parser


@pytest.fixture
def maf(tmp_path):
    path = tmp_path / "sample.maf"
    path.write_text("Chromosome\tStart_Position\n")
    return path


class TestCollect:
    def test_file_arguments_become_paths(self, maf):
        args = build_parser().parse_args(["-i", str(maf)])

        assert collect_run_parameters(args)["infile"] == [str(maf)]

    def test_stdout_is_not_recorded(self, maf):
        """There is no path to write down, and "<stdout>" is not one."""
        args = build_parser().parse_args(["-i", str(maf)])

        assert "outfile" not in collect_run_parameters(args)

    def test_named_output_file_is_recorded(self, tmp_path, maf):
        out = tmp_path / "out.tsv"
        args = build_parser().parse_args(["-i", str(maf), "-o", str(out)])

        assert collect_run_parameters(args)["outfile"] == str(out)

    def test_plain_values_pass_through(self, maf):
        args = build_parser().parse_args(["-i", str(maf), "-g", "hg38", "-n", "42"])
        parameters = collect_run_parameters(args)

        assert parameters["genome"] == "hg38"
        assert parameters["nsamples"] == 42

    def test_the_parameter_flags_are_not_themselves_parameters(self, tmp_path, maf):
        """Replaying --params-out would overwrite the record it came from."""
        args = build_parser().parse_args(
            ["-i", str(maf), "--params-out", str(tmp_path / "r.json"), "--params-in", "x.json"]
        )
        parameters = collect_run_parameters(args)

        assert "params_out" not in parameters
        assert "params_in" not in parameters

    def test_every_value_is_json_serializable(self, maf):
        args = build_parser().parse_args(["-i", str(maf)])

        json.dumps(collect_run_parameters(args))  # must not raise


class TestWriteAndLoad:
    def test_round_trip(self, tmp_path, maf):
        args = build_parser().parse_args(["-i", str(maf), "-g", "hg38"])
        path = tmp_path / "run.json"

        write_run_parameters(args, "profile", path)
        command, arguments = load_run_parameters(path)

        assert command == "profile"
        assert arguments["genome"] == "hg38"
        assert arguments["infile"] == [str(maf)]

    def test_unwritable_destination_warns_rather_than_raising(self, tmp_path, maf):
        args = build_parser().parse_args(["-i", str(maf)])

        assert write_run_parameters(args, "profile", tmp_path / "no" / "such" / "r.json") is None

    def test_invalid_json_is_reported_clearly(self, tmp_path):
        path = tmp_path / "bad.json"
        path.write_text("not json")

        with pytest.raises(ValueError, match="not valid JSON"):
            load_run_parameters(path)

    def test_a_foreign_json_file_is_rejected(self, tmp_path):
        path = tmp_path / "other.json"
        path.write_text('{"something": "else"}')

        with pytest.raises(ValueError, match="does not look like"):
            load_run_parameters(path)

    def test_missing_file_is_reported_clearly(self, tmp_path):
        with pytest.raises(ValueError, match="Could not read"):
            load_run_parameters(tmp_path / "absent.json")

    def test_version_mismatch_warns_but_still_loads(self, tmp_path, caplog):
        path = tmp_path / "old.json"
        path.write_text(
            json.dumps(
                {"mutagene_version": "0.0.1", "command": "profile", "arguments": {"genome": "hg19"}}
            )
        )

        command, arguments = load_run_parameters(path)

        assert (command, arguments["genome"]) == ("profile", "hg19")
        assert "0.0.1" in caplog.text


class TestReplay:
    """The command line must beat the file, not the other way round."""

    def replay(self, arguments, argv):
        parser = build_parser()
        parser.set_defaults(**arguments)
        args = parser.parse_args(argv)
        convert_recorded_lists(parser, args)
        return args

    def test_recorded_values_are_used_when_no_argument_is_given(self, maf):
        args = self.replay({"genome": "hg38", "infile": [str(maf)]}, [])

        assert args.genome == "hg38"

    def test_command_line_overrides_the_file(self, maf):
        args = self.replay({"genome": "hg38", "infile": [str(maf)]}, ["-g", "mm10"])

        assert args.genome == "mm10"

    def test_recorded_paths_are_reopened_as_files(self, maf):
        args = self.replay({"infile": [str(maf)]}, [])

        assert not isinstance(args.infile[0], str), "a path was left unopened"
        assert args.infile[0].read().startswith("Chromosome")

    def test_command_line_files_are_not_reopened_from_the_record(self, tmp_path, maf):
        other = tmp_path / "other.maf"
        other.write_text("OTHER\n")

        args = self.replay({"infile": [str(maf)]}, ["-i", str(other)])

        assert [f.name for f in args.infile] == [str(other)]

    def test_typed_scalars_are_converted(self, maf):
        args = self.replay({"nsamples": "42", "infile": [str(maf)]}, [])

        assert args.nsamples == 42

    def test_a_recorded_output_file_is_not_opened_when_overridden(self, tmp_path, maf):
        """Opening it eagerly would truncate a previous run's output."""
        previous = tmp_path / "previous.tsv"
        previous.write_text("EARLIER RESULTS")
        new = tmp_path / "new.tsv"

        self.replay({"outfile": str(previous), "infile": [str(maf)]}, ["-o", str(new)])

        assert previous.read_text() == "EARLIER RESULTS"


class TestPeekParamsIn:
    def test_finds_the_flag_after_the_subcommand(self):
        assert peek_params_in(["profile", "--params-in", "r.json", "-g", "hg19"]) == "r.json"

    def test_accepts_the_equals_form(self):
        assert peek_params_in(["profile", "--params-in=r.json"]) == "r.json"

    def test_absent_flag_gives_none(self):
        assert peek_params_in(["profile", "-g", "hg19"]) is None

    def test_unknown_arguments_do_not_upset_it(self):
        assert peek_params_in(["--nonsense", "--params-in", "r.json"]) == "r.json"
