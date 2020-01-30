import argparse
import signal
import sys

import mutagene
from mutagene.version import __version__

# from mutagene.io.profile import format_profile
# from mutagene.io.mutations_profile import read_VCF_profile
# from mutagene.io.decomposition import write_decomposition

from mutagene.cli.fetch_menu import FetchMenu
from mutagene.cli.profile_menu import ProfileMenu
from mutagene.cli.motif_menu import MotifMenu
from mutagene.cli.signature_menu import SignatureMenu
from mutagene.cli.rank_menu import RankMenu
from mutagene.cli.benchmark_menu import BenchmarkMenu

import logging

logger = logging.getLogger(__name__)

genome_error_message = 'requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more'


class MutaGeneApp(object):
    def __init__(self):
        signal.signal(signal.SIGINT, self.signal_handler)
        # ignore BrokenPipeError: [Errno 32] Broken pipe which occurs when using less or head
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

        parser = argparse.ArgumentParser(
            prog="mutagene",
            description='MutaGene version {} - Analysis of mutational processes and driver mutations'.format(mutagene.__version__),
            formatter_class=argparse.RawDescriptionHelpFormatter,
            add_help=False,
            # usage="%(prog)s [options]",
            # formatter_class=argparse.RawTextHelpFormatter
        )
        # parser._positionals.title = 'Positional arguments'
        parser._optionals.title = 'Global optional arguments'
        parser.add_argument('-v', '--verbose', action='count', default=0, help='Print additional messages (-v, -vv)')
        parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__, help='Show version and exit')
        parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')

        subparsers = parser.add_subparsers(
            help="",
            metavar='{fetch, profile, rank, motif, signature}',
            description="""\
        fetch - Load data such as genomes and cancer datasets from demote sources (alias: download)
        profile - Create a mutational profile given a sample with mutations
        rank - Predict driver mutations by ranking observed mutations with respect to their expected mutability
        motif - Test samples for presence of mutational motifs
        signature - Identify activity of existing mutational signatures in samples or derive new signatures (aliases: identify, decompose)\
            """, dest='command', title='Choose MutaGene subpackage')

        parser_mapping = {
            'fetch': {'class': FetchMenu, 'aliases': ['download']},
            'profile': {'class': ProfileMenu, 'aliases': []},
            'rank': {'class': RankMenu, 'aliases': ['driver']},
            'motif': {'class': MotifMenu, 'aliases': []},
            'signature': {'class': SignatureMenu, 'aliases': ['identify', 'decompose']},
            'benchmark': {'class': BenchmarkMenu, 'aliases': []},
        }

        for command, menu in parser_mapping.items():
            # initialize parser object for each subparser
            parser_mapping[command]['parser'] = menu['class'](
                subparsers.add_parser(
                    command,
                    add_help=True,
                    aliases=menu['aliases'],
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                ))

        args = parser.parse_args()

        if not args.command:
            parser.print_help()
            sys.exit(1)

        # LOGGER INIT
        levels = [logging.WARNING, logging.INFO, logging.DEBUG]
        level = levels[min(len(levels) - 1, args.verbose)]  # capped to number of levels
        logging.basicConfig(level=level,
                            format="%(levelname)s %(message)s")

        # Calling callback method for parser object by command name or alias
        parser_class = None
        if args.command in parser_mapping:
            parser_class = parser_mapping[args.command]['parser']
        else:
            for command, mapping in parser_mapping.items():
                if args.command in mapping['aliases']:
                    parser_class = mapping['parser']
                    break

        # should not happen if we have a correct name or an alias for a command
        if not parser_class:
            parser.print_help()
            sys.exit(1)

        parser_class.callback(args)

    @classmethod
    def signal_handler(cls, signal, frame):
        # logger.warning('Interrupted')
        sys.exit(0)


if __name__ == '__main__':
    MutaGeneApp()
