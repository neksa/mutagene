import argparse
import sys
import logging

from mutagene.io.context_window import read_MAF_with_context_window
from mutagene.motifs import identify_motifs
from mutagene.io.motifs import write_motif_matches
from mutagene.motifs import motifs as list_of_motifs


logger = logging.getLogger(__name__)
genome_error_message = """requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more
                        Use mutagene fetch to download genome assemblies"""


class MotifMenu(object):
    def __init__(self, parser):

        parser.description = "Motif function requires: mutagene motif <action (search or list)>, if search is specified, infile & genome are also required"
        parser.epilog = """
                        Examples:
                        # search for the presence of the C[A>T] motif in sample1.maf using hg19
                        mutagene motif search -i sample1.maf -g hg19 -m 'C[A>T]'
                        # search in sample2.vcf for all preidentified motifs in mutagene using hg18
                        mutagene motif search -i sample2.vcf -g hg18
                        """

        parser.add_argument('action', choices=['search', 'list'], help="search for a motif, list all predefined motifs")

        parser.add_argument("--infile", "-i", help="Input file in MAF or VCF format with one or multiple samples",
                            type=argparse.FileType('r'))
        parser.add_argument('--genome', "-g", help="Location of genome assembly file in 2bit format", type=str)

        parser.add_argument("--motif", "-m",
                            help="Motif to search for, use the 'R[C>T]GY' syntax for the motif. Use quotes", type=str)
        parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                            help="Name of output file, will be generated in TSV format")
        parser.add_argument('--window-size', "-w", help="Context window size for motif search, default setting is 50",
                            type=int, default=50)
        parser.add_argument('--strand', "-s",
                            help="Transcribed strand (+), non-transcribed (-), any (=), or all (+-= default) ",
                            type=str, default='+-=', choices=['+', '-', '=', '+-='])

    @classmethod
    def search(cls, args):
        if not args.infile:
            logger.warning("Provide input file in VCF or MAF format (-i) and a corresponding genome assembly (-g)")
            return
        if not args.genome:
            logger.warning(genome_error_message)
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
                    "Mutagene motif search failed because motif cannot be processed. Check to make sure motif input is in quotes")
                return

            logger.info("Searching for motif {}".format(custom_motif))

        if args.window_size > 250 or args.window_size < 1:
            logger.warning('window-size should be between 1 and 250 nucleotides')
            return

        mutations, mutations_with_context, processing_stats = read_MAF_with_context_window(args.infile, args.genome,
                                                                                           args.window_size)
        if len(mutations_with_context) == 0:
            logger.warning("No mutations loaded")

        matching_motifs = identify_motifs(mutations_with_context, custom_motif,
                                          args.strand) if mutations_with_context is not None else []

        if len(matching_motifs) == 0:
            logger.warning("No significant motif matches found")
        else:
            write_motif_matches(args.outfile, matching_motifs)

    @classmethod
    def list(cls, args):

        if args.genome:
            logger.warning("Genome argument not accepted for motif list")

        if args.infile:
            logger.warning("Infile argument not accepted for motif list")

        if args.motif:
            logger.warning("Motif argument not accepted for motif list")

        if args.window_size != 50:
            logger.warning("Window size argument not accepted for motif list")

        if args.strand != '*':
            logger.warning("Strand argument not accepted for motif list")

        if args.outfile != sys.stdout:
            logger.warning("Outfile: argument not accepted for motif list")

        for m in list_of_motifs:
            print("{:20}\t{}".format(m['name'], m['logo']))

    @classmethod
    def callback(cls, args):
        # print('MotifMenu', args.action)
        getattr(cls, args.action)(args)
