import argparse
import sys
import logging

from mutagene.io.context_window import read_MAF_with_context_window, read_VCF_with_context_window
from mutagene.motifs import identify_motifs
from mutagene.io.motifs import write_motif_matches, get_known_motifs


logger = logging.getLogger(__name__)
genome_error_message = """requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more
                        Use mutagene fetch to download genome assemblies"""


class MotifMenu(object):
    def __init__(self, parser):
        parser.description = ""
        parser.epilog = """
Examples:
# search in sample2.vcf for all preidentified motifs in mutagene using hg19
mutagene motif --infile sample2.vcf --input-format VCF --genome hg19

# search for the presence of the C[A>T] motif in sample1.maf using hg19 not checking for strand-specificity
mutagene motif --infile sample1.maf --input-format MAF --genome hg19 --motif 'C[A>T]' --strand A
        """

        ###################################################################
        required_group = parser.add_argument_group('Required arguments')
        required_group.add_argument(
            "--infile", "-i", help="Input file in MAF or VCF format with one or multiple samples",
            type=argparse.FileType('r'))
        required_group.add_argument(
            '--genome', "-g", help="Location of genome assembly file in 2bit format", type=str)

        ###################################################################
        optional_group = parser.add_argument_group('Optional arguments')
        optional_group.add_argument('--input-format', "-f", help="Input format: MAF, VCF", type=str, choices=['MAF', 'VCF'], default='MAF')
        optional_group.add_argument(
            "--motif", "-m",
            help="Motif to search for, use the 'R[C>T]GY' syntax for the motif. Use quotes", type=str)
        optional_group.add_argument(
            '--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout,
            help="Name of output file, will be generated in TSV format")

        # for backwards compatibility with 0.8.X add a hidden action that would just take anything as a valid input
        optional_group.add_argument('action', nargs="?", metavar="")

        ###################################################################
        advanced_group = parser.add_argument_group('Advanced arguments')
        advanced_group.add_argument(
            '--window-size', "-w", help="Context window size for motif search, default setting is 50",
            type=int, default=50)
        advanced_group.add_argument(
            '--strand', "-s",
            help="Transcribed strand (T), non-transcribed (N), any (A), or all (TNA default) ",
            type=str, default='TNA', choices=['T', 'N', 'A', 'TNA'])
        advanced_group.add_argument(
            '--threshold', "-t",
            help="Significance threshold for qvalues, default value=0.05",
            type=float, default=0.05)
        advanced_group.add_argument(
            '--save-motif-matches',
            help="Save mutations in matching motifs to a BED file",
            type=argparse.FileType('w'), default=None)
        advanced_group.add_argument(
            '--test',
            help="Statistical test to use",
            type=str, default='Fisher', choices=['Fisher', 'Chi2'])

        self.parser = parser

    @classmethod
    def search(cls, args):
        if not args.infile:
            logger.warning("Provide input file in VCF or MAF format (-i) and a corresponding genome assembly (-g)")
            return
        if not args.genome:
            logger.warning(genome_error_message)
            return
        if args.threshold > 1.0 or args.threshold < 0.0:
            logger.warning("The threshold value should be between 0.0 and 1.0, inclusive")
            return
        if not args.motif:
            logger.info("Searching for predefined motifs")
            logger.info("Will reduce statistical power because of multiple comparisons")
            custom_motif = None
        else:
            custom_motif = args.motif
            custom_motif = custom_motif.replace('-', '>')
            custom_motif = custom_motif.replace('/', '>')
            custom_motif = custom_motif.replace('.', '>')
            custom_motif = custom_motif.replace('->', '>')

            if ">" not in custom_motif or "]" not in custom_motif:
                logger.warning(
                    "Mutagene motif search failed because motif cannot be processed."
                    "Check to make sure motif input is in quotes")
                return

            logger.info("Searching for motif {}".format(custom_motif))

        if args.window_size > 5000 or args.window_size < 1:
            logger.warning('window-size should be between 1 and 5000 nucleotides')
            return

        try:
            if args.input_format == 'VCF':
                mutations, mutations_with_context, processing_stats = read_VCF_with_context_window(
                    args.infile, args.genome, args.window_size)
            elif args.input_format == 'MAF':
                mutations, mutations_with_context, processing_stats = read_MAF_with_context_window(
                    args.infile, args.genome, args.window_size)
        except ValueError as e:
            logger.warning('Not able to parse input file in {} format: {}. You can specify a different format with --input-format (-f)'.format(args.input_format, e))
            sys.exit(1)

        if len(mutations_with_context) == 0:
            logger.warning("No mutations loaded")

        #######################
        # Performance PROFILING
        # import cProfile
        # import pstats
        # pr = cProfile.Profile()
        # pr.enable()

        matching_motifs = identify_motifs(
            samples_mutations=mutations_with_context,
            custom_motif=custom_motif,
            strand=args.strand,
            threshold=args.threshold,
            dump_matches=args.save_motif_matches,
            stat_type=args.test) if mutations_with_context is not None else []

        #######################
        # Performance PROFILING
        # pr.disable()
        # p = pstats.Stats(pr)
        # p.sort_stats('ncalls').print_stats(20)  # strip_dirs()

        if len(matching_motifs) == 0:
            logger.warning("No significant motif matches found")
        else:
            write_motif_matches(args.outfile, matching_motifs)

    def list(self, args):
        """ Prints the list of known motifs bundled with the package"""
        print("\nThe list of mutational motifs that will be tested by default:")
        for m in get_known_motifs():
            print("{:20}\t{}".format(m['name'], m['logo']))
        print("any custom motif can be specified with --motif (-m)")

    def callback(self, args):
        # getattr(cls, args.action)(args)
        if args.infile:
            MotifMenu.search(args)
        else:
            self.parser.print_usage()
            self.list(args)
            sys.exit(1)
