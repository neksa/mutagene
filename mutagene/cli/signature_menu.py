import argparse
import sys
import logging

from mutagene.io.profile import read_signatures
from mutagene.profiles.profile import get_mutational_profile
from mutagene.io.mutations_profile import read_auto_profile
from mutagene.signatures.identify import decompose_mutational_profile_counts
from mutagene.io.decomposition import write_decomposition


logger = logging.getLogger(__name__)
genome_error_message = 'requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more'


class SignatureMenu(object):
    def __init__(self, parser):
        parser.add_argument('action', choices=['identify', ])  # 'new'])
        parser.add_argument("--signatures", "-s", choices=[5, 10, 30], help="Collection of signatures to use", type=int)
        parser.add_argument("--infile", "-i", help="Input file in VCF or MAF format", type=argparse.FileType('r'))
        parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
        parser.add_argument('--genome', "-g", help="Location of genome assembly file in 2bit format", type=str)

    @classmethod
    def identify(cls, args):
        if not args.infile:
            logger.warning("Provide input file in VCF or MAF format (-i) and a corresponding genome assembly (-g)")
            return
        if not args.genome:
            logger.warning(genome_error_message)
            return
        if not args.signatures:
            logger.warning("Set of signatures required. Use 5 and 10 for MUTAGENE-5 and MUTAGENE-10. Use 30 for COSMIC-30")
            return

        # mutations, processing_stats = read_VCF_profile(args.infile, asm=args.genome)
        mutations, processing_stats = read_auto_profile(args.infile, fmt='auto', asm=args.genome)
        W, signature_names = read_signatures(int(args.signatures))
        # print(W, signature_names)
        profile = get_mutational_profile(mutations, counts=True)

        _, _, results = decompose_mutational_profile_counts(
            profile,
            (W, signature_names),
            'MLEZ',
            debug=False,
            others_threshold=0.0)
        write_decomposition(args.outfile, results, signature_names)

    @classmethod
    def callback(cls, args):
        # print('SignatureMenu', args.action)
        getattr(cls, args.action)(args)
