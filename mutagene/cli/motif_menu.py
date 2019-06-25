import argparse
import sys
import logging

from mutagene.io.context_window import read_MAF_with_context_window
from mutagene.motifs.motifs import identify_motifs
from mutagene.io.motifs import write_motif_matches


logger = logging.getLogger(__name__)
genome_error_message = 'requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more'


class MotifMenu(object):
    def __init__(self, parser):
        parser.add_argument('action', choices=['search'])
        parser.add_argument("--motif", "-m", help="Motif to search for, use the 'R[C>T]GY' syntax for the motif. Use quotes", type=str)
        parser.add_argument("--infile", "-i", help="Input file in MAF or VCF format with one or multiple samples", type=argparse.FileType('r'))
        parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
        parser.add_argument('--genome', "-g", help="Location of genome assembly file in 2bit format", type=str)
        parser.add_argument('--window-size', "-w", help="Context window size for motif search, default setting is 50", type=int, default=50)
        parser.add_argument('--strand', "-s", help="Transcribed strand (+), non-transcribed (-), or both (*): the default setting", type=str, default='*', choices=['*', '+', '-'])

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
            custom_motif = None
        else:
            custom_motif = args.motif
            custom_motif = custom_motif.replace('-', '>')
            custom_motif = custom_motif.replace('/', '>')
            custom_motif = custom_motif.replace('.', '>')
            custom_motif = custom_motif.replace('->', '>')
            logger.info("Searching for motif {}".format(custom_motif))

        if args.window_size > 250 or args.window_size < 1:
            logger.warning('window-size should be between 1 and 250 nucleotides')
            return

        mutations, mutations_with_context, processing_stats = read_MAF_with_context_window(args.infile, args.genome, args.window_size)
        matching_motifs = identify_motifs(mutations_with_context, custom_motif, args.strand) if mutations_with_context is not None else []
        if len(matching_motifs) == 0:
            logger.warning("No significant motif matches found")
        else:
            write_motif_matches(args.outfile, matching_motifs)

    @classmethod
    def callback(cls, args):
        # print('MotifMenu', args.action)
        getattr(cls, args.action)(args)
