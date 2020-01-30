import argparse
import sys
import logging
from mutagene.profiles.profile import calc_profile


logger = logging.getLogger(__name__)
genome_error_message = """requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more
                          Use mutagene fetch to download genome assemblies"""


class ProfileMenu(object):
    def __init__(self, parser):
        parser.add_argument("--infile", "-i", nargs='*', help="Input file format", type=argparse.FileType('r'))
        parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                            help="Name of output file, will be generated in TSV format")
        parser.add_argument('--genome', "-g", help="Location of genome assembly file", type=str)
        parser.add_argument('--input-format', "-f", help="Input format: auto, MAF, VCF", type=str, default='auto')

        # for backwards compatibility with 0.8.X add a hidden action that would just take anything as a valid input
        parser.add_argument('action', nargs="?", metavar="")

    def callback(self, args):
        # print('ProfileMenu', args.action)
        self.calculate(args)

    def calculate(self, args):
        # print("Calculating...")
        if not args.infile:
            logger.warning("Provide input file in VCF or MAF format (-i) and a corresponding genome assembly (-g)")
            return
        if not args.genome:
            logger.warning(genome_error_message)
            return
        calc_profile(args.infile, args.outfile, args.genome, args.input_format)
