"""Recording a run's parameters and replaying them (GitHub issue #45).

``--params-out FILE`` writes what a run was actually asked to do; ``--params-in
FILE`` reads it back as defaults, so a result can be reproduced from the record
rather than from shell history. Both are opt-in: nothing is written unless
asked.

Replay works through ``parser.set_defaults`` rather than by overwriting the
parsed namespace. That is what makes an argument given on the command line beat
the same argument in the file: argparse only falls back to a default when the
argument is absent, and it applies the argument's own ``type`` to a string
default, so a recorded path is reopened as a file exactly as it would have been
from the command line.
"""

import argparse
import json
import logging
from pathlib import Path

from mutagene.version import __version__

logger = logging.getLogger(__name__)

# None of these describe what the run did. ``params_out`` would make a replay
# overwrite the record it came from; ``command`` is stored separately and is
# never a default; ``action`` is the hidden 0.8.x compatibility positional.
_EXCLUDED = frozenset({"params_in", "params_out", "command", "action"})

# Distinguishes "this value cannot be written down" from a genuine None.
_UNRECORDABLE = object()


def add_run_parameter_arguments(parser):
    """Attach --params-out/--params-in to a subcommand parser."""
    group = parser.add_argument_group("Reproducibility arguments")
    group.add_argument(
        "--params-out",
        metavar="FILE",
        type=str,
        help="Write the parameters of this run to FILE in JSON format",
    )
    group.add_argument(
        "--params-in",
        metavar="FILE",
        type=str,
        help="Take parameters from a FILE written by --params-out. "
        "Arguments given on the command line override it",
    )


def _describe(value):
    """Reduce an argument value to something JSON holds and argparse can re-read."""
    if isinstance(value, bool | int | float | str) or value is None:
        return value

    if isinstance(value, list | tuple):
        # nargs='*' arguments, such as `profile --infile`, arrive as a list of
        # open files. Each element has to be reduced on its own; stringifying
        # the list gives a repr that no replay can reopen.
        described = [_describe(item) for item in value]
        return _UNRECORDABLE if any(d is _UNRECORDABLE for d in described) else described

    name = getattr(value, "name", None)
    if name is not None:
        # An open file from argparse.FileType. Record where it points rather
        # than the handle: a path is what makes the record replayable. Streams
        # such as stdout report a name like "<stdout>" and have no path to keep.
        return _UNRECORDABLE if str(name).startswith("<") else str(name)

    return str(value)


def collect_run_parameters(args):
    """The run's arguments as plain JSON-able values, in a stable order."""
    parameters = {}
    for key, value in sorted(vars(args).items()):
        if key in _EXCLUDED:
            continue
        described = _describe(value)
        if described is _UNRECORDABLE:
            continue
        parameters[key] = described
    return parameters


def write_run_parameters(args, command, path):
    """Write the run's parameters to *path*; returns the path, or None on failure."""
    record = {
        "mutagene_version": __version__,
        "command": command,
        "arguments": collect_run_parameters(args),
    }
    try:
        Path(path).write_text(json.dumps(record, indent=2, sort_keys=False) + "\n")
    except OSError as e:
        logger.warning(f"Could not write run parameters to {path}: {e}")
        return None
    logger.info(f"Run parameters written to {path}")
    return path


def load_run_parameters(path):
    """Read a run-parameter file; returns ``(command, arguments)``.

    Raises ValueError with a usable message when the file is not one of ours.
    """
    try:
        record = json.loads(Path(path).read_text())
    except OSError as e:
        raise ValueError(f"Could not read run parameters from {path}: {e}") from None
    except json.JSONDecodeError as e:
        raise ValueError(f"{path} is not valid JSON: {e}") from None

    if not isinstance(record, dict) or "arguments" not in record:
        raise ValueError(f"{path} does not look like a mutagene run-parameter file")

    arguments = record["arguments"]
    if not isinstance(arguments, dict):
        raise ValueError(f"'arguments' in {path} must be an object")

    version = record.get("mutagene_version")
    if version and version != __version__:
        logger.warning(
            f"{path} was written by mutagene {version}, running {__version__}. "
            "Arguments may have changed meaning since"
        )

    return record.get("command"), arguments


def convert_recorded_lists(subparser, args):
    """Apply each argument's own ``type`` to list values that came from a file.

    argparse converts a *string* default with the argument's type when the
    argument is absent, which is what makes a recorded path reopen as a file.
    It does not do this element by element for a list default, so a recorded
    ``nargs='*'`` argument such as ``profile --infile`` would arrive as a list
    of plain strings and be read a character at a time.

    A list still holding strings can only have come from the defaults: had the
    argument been given on the command line, argparse would already have
    converted it. So there is nothing to detect here, and files are still not
    opened until the value is actually needed.
    """
    for action in subparser._actions:
        if action.type is None or action.type is str or action.dest is None:
            continue
        value = getattr(args, action.dest, None)
        if not isinstance(value, list):
            continue
        setattr(
            args,
            action.dest,
            [action.type(item) if isinstance(item, str) else item for item in value],
        )


def peek_params_in(argv):
    """Find --params-in in *argv* before the real parser is built.

    The file supplies defaults, and defaults have to be in place before parsing
    for the command line to be able to override them.
    """
    bootstrap = argparse.ArgumentParser(add_help=False)
    bootstrap.add_argument("--params-in")
    try:
        known, _ = bootstrap.parse_known_args(argv)
    except SystemExit:
        # Malformed argv; let the real parser produce the error message.
        return None
    return known.params_in
