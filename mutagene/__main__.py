import argparse
import logging
import signal
import sys

import mutagene
from mutagene.cli.benchmark_menu import BenchmarkMenu
from mutagene.cli.demo_menu import DemoMenu

# from mutagene.io.profile import format_profile
# from mutagene.io.mutations_profile import read_VCF_profile
# from mutagene.io.decomposition import write_decomposition
from mutagene.cli.fetch_menu import FetchMenu
from mutagene.cli.motif_menu import MotifMenu
from mutagene.cli.profile_menu import ProfileMenu
from mutagene.cli.rank_menu import RankMenu
from mutagene.cli.run_params import (
    add_run_parameter_arguments,
    convert_recorded_lists,
    load_run_parameters,
    peek_params_in,
    write_run_parameters,
)
from mutagene.cli.serve_menu import ServeMenu  # lazy import in serve_menu.callback
from mutagene.cli.signature_menu import SignatureMenu
from mutagene.version import __version__

logger = logging.getLogger(__name__)

genome_error_message = (
    "requires genome name argument -g hg19, hg38, mm10, "
    "see http://hgdownload.soe.ucsc.edu/downloads.html for more. "
    "Use mutagene fetch to download genome assemblies"
)


class MutaGeneApp:
    def __init__(self):
        signal.signal(signal.SIGINT, self.signal_handler)
        # ignore BrokenPipeError: [Errno 32] Broken pipe which occurs when using less or head
        if hasattr(signal, "SIGPIPE"):
            signal.signal(signal.SIGPIPE, signal.SIG_DFL)

        parser = argparse.ArgumentParser(
            prog="mutagene",
            description=f"MutaGene version {mutagene.__version__} - Analysis of mutational processes and driver mutations",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            add_help=False,
            # usage="%(prog)s [options]",
            # formatter_class=argparse.RawTextHelpFormatter
        )
        # parser._positionals.title = 'Positional arguments'
        parser._optionals.title = "Global optional arguments"
        parser.add_argument(
            "-v", "--verbose", action="count", default=0, help="Print additional messages (-v, -vv)"
        )
        parser.add_argument(
            "-V",
            "--version",
            action="version",
            version="%(prog)s " + __version__,
            help="Show version and exit",
        )
        parser.add_argument(
            "-h",
            "--help",
            action="help",
            default=argparse.SUPPRESS,
            help="Show this help message and exit",
        )

        subparsers = parser.add_subparsers(
            dest="command",
            title="Commands",
            metavar="COMMAND",
        )

        parser_mapping = {
            "fetch": {
                "help": "Download genomes and reference data",
                "class": FetchMenu,
                "aliases": ["download"],
            },
            "profile": {
                "help": "Count mutations by sequence context",
                "class": ProfileMenu,
                "aliases": [],
            },
            "rank": {
                "help": "Rank mutations by expected mutability",
                "class": RankMenu,
                "aliases": ["driver"],
            },
            "motif": {"help": "Test for motif enrichment", "class": MotifMenu, "aliases": []},
            "signature": {
                "help": "Fit known mutational signatures",
                "class": SignatureMenu,
                "aliases": ["identify", "decompose"],
            },
            "benchmark": {
                "help": "Benchmark fits on simulated profiles",
                "class": BenchmarkMenu,
                "aliases": [],
            },
            "serve": {"help": "Start the local web interface", "class": ServeMenu, "aliases": []},
            "demo": {"help": "Run the offline synthetic example", "class": DemoMenu, "aliases": []},
        }

        for command, menu in parser_mapping.items():
            # initialize parser object for each subparser
            subparser = subparsers.add_parser(
                command,
                help=menu["help"],
                aliases=menu["aliases"],
                formatter_class=argparse.RawDescriptionHelpFormatter,
            )
            add_run_parameter_arguments(subparser)
            parser_mapping[command]["subparser"] = subparser
            parser_mapping[command]["parser"] = menu["class"](subparser)

        # --params-in has to be resolved before parsing, because it supplies
        # defaults and argparse only falls back to a default when the argument
        # is absent. Applying it afterwards would let the file override the
        # command line, which is the wrong way round.
        recorded_command = None
        params_in = peek_params_in(sys.argv[1:])
        if params_in:
            try:
                recorded_command, recorded_arguments = load_run_parameters(params_in)
            except ValueError as e:
                logger.error(str(e))
                sys.exit(1)
            if recorded_command not in parser_mapping:
                logger.error(
                    f"{params_in} records the command {recorded_command!r}, "
                    "which this version of mutagene does not have"
                )
                sys.exit(1)
            parser_mapping[recorded_command]["subparser"].set_defaults(**recorded_arguments)

        args = parser.parse_args()

        if not args.command:
            parser.print_help()
            sys.exit(1)

        # LOGGER INIT
        levels = [logging.WARNING, logging.INFO, logging.DEBUG]
        level = levels[min(len(levels) - 1, args.verbose)]  # capped to number of levels
        logging.basicConfig(level=level, format="%(levelname)s %(message)s")

        # Calling callback method for parser object by command name or alias
        parser_class = None
        canonical_command = None
        if args.command in parser_mapping:
            parser_class = parser_mapping[args.command]["parser"]
            canonical_command = args.command
        else:
            for command, mapping in parser_mapping.items():
                if args.command in mapping["aliases"]:
                    parser_class = mapping["parser"]
                    canonical_command = command
                    break

        # should not happen if we have a correct name or an alias for a command
        if not parser_class:
            parser.print_help()
            sys.exit(1)

        if recorded_command and recorded_command != canonical_command:
            logger.error(
                f"{params_in} holds parameters for '{recorded_command}', "
                f"not '{canonical_command}'"
            )
            sys.exit(1)

        if recorded_command:
            convert_recorded_lists(parser_mapping[canonical_command]["subparser"], args)

        try:
            parser_class.callback(args)
        finally:
            # Written in a finally so a run that fails still leaves a record of
            # what it was asked to do, and written last so it captures values
            # the callback resolved, such as the full path of the genome.
            if args.params_out:
                write_run_parameters(args, canonical_command, args.params_out)

    @classmethod
    def signal_handler(cls, signal, frame):
        sys.exit(130)


def main():
    MutaGeneApp()


if __name__ == "__main__":
    main()
