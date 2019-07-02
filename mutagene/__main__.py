import argparse
import signal
import sys

import mutagene
from mutagene.version import __version__

# from mutagene.io.profile import format_profile
# from mutagene.io.mutations_profile import read_VCF_profile
# from mutagene.io.decomposition import write_decomposition

from mutagene.cli.fetch_menu import FetchMenu
from mutagene.cli.motif_menu import MotifMenu
from mutagene.cli.profile_menu import ProfileMenu
from mutagene.cli.signature_menu import SignatureMenu
from mutagene.cli.rank_menu import RankMenu

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
            # usage="%(prog)s [options]",
            # formatter_class=argparse.RawTextHelpFormatter
        )
        parser.add_argument('-v', '--verbose', action='count', default=0)
        parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

        subparsers = parser.add_subparsers(help='sub-command help', dest='command')

        parser_mapping = {
            'fetch': FetchMenu,
            'profile': ProfileMenu,
            'rank': RankMenu,
            'motif': MotifMenu,
            'signature': SignatureMenu,
        }

        for command, menu in parser_mapping.items():
            menu(subparsers.add_parser(command))

        args = parser.parse_args()

        if not args.command:
            parser.print_help()
            sys.exit(1)

        # LOGGER INIT
        levels = [logging.WARNING, logging.INFO, logging.DEBUG]
        level = levels[min(len(levels) - 1, args.verbose)]  # capped to number of levels
        logging.basicConfig(level=level,
                            format="%(levelname)s %(message)s")

        parser_mapping[args.command].callback(args)

    @classmethod
    def signal_handler(cls, signal, frame):
        # logger.warning('Interrupted')
        sys.exit(0)


if __name__ == '__main__':
    MutaGeneApp()
