import argparse
import logging
import sys

from mutagene.cli.genome_resolver import GenomeNotFoundError, resolve_genome
from mutagene.io.variant_filter import add_filter_argument
from mutagene.profiles.profile import calc_profile

logger = logging.getLogger(__name__)
genome_error_message = """requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.soe.ucsc.edu/downloads.html for more
                          Use mutagene fetch to download genome assemblies"""


class ProfileMenu:
    def __init__(self, parser):
        parser.add_argument(
            "--infile", "-i", nargs="*", help="Input file format", type=argparse.FileType("r")
        )
        parser.add_argument(
            "--outfile",
            "-o",
            nargs="?",
            type=argparse.FileType("w"),
            default=sys.stdout,
            help="Name of output file, will be generated in TSV format",
        )
        parser.add_argument("--genome", "-g", help="Location of genome assembly file", type=str)
        parser.add_argument(
            "--input-format", "-f", help="Input format: auto, MAF, VCF", type=str, default="auto"
        )

        # for backwards compatibility with 0.8.X add a hidden action that would just take anything as a valid input
        add_filter_argument(parser)

        parser.add_argument("action", nargs="?", metavar="")

    def callback(self, args):
        # print('ProfileMenu', args.action)
        self.calculate(args)

    def calculate(self, args):
        # print("Calculating...")
        if not args.infile:
            logger.warning(
                "Provide input file in VCF or MAF format (-i) and a corresponding genome assembly (-g)"
            )
            return
        if not args.genome:
            logger.warning(genome_error_message)
            return
        try:
            args.genome = resolve_genome(args.genome)
        except GenomeNotFoundError as e:
            logger.error(str(e))
            sys.exit(1)

        try:
            calc_profile(
                args.infile,
                args.outfile,
                args.genome,
                args.input_format,
                keep_filtered=args.keep_filtered,
            )
        except (ValueError, OSError) as e:
            # An unreadable input is the user's problem to fix, not a defect to
            # report, so say what is wrong instead of printing a traceback.
            # rank and signature already did this; profile did not.
            logger.error(
                f"Could not read the input in {args.input_format} format: {e}\n"
                "Check the file, or name a different format with --input-format (-f)"
            )
            if logger.root.level == logging.DEBUG:
                raise
            sys.exit(1)
