import argparse
import sys
import logging

from mutagene.io.profile import read_signatures
from mutagene.profiles.profile import get_multisample_mutational_profile, get_mutational_profile
from mutagene.profiles.profile import generate_resampled_profiles
from mutagene.io.mutations_profile import read_auto_profile
from mutagene.io.context_window import read_MAF_with_context_window
from mutagene.io.context_window import read_TCGI_with_context_window
from mutagene.signatures.identify import decompose_mutational_profile_counts
from mutagene.signatures.identify import IDENTIFY_MIN_FUNCTIONS
from mutagene.io.decomposition import write_multisample_decomposition
from mutagene.io.decomposition import write_decomposition
from mutagene.io.decomposition import write_bootstrap_decomposition


logger = logging.getLogger(__name__)

genome_error_message = """requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more
                          Use mutagene fetch to download genome assemblies"""


class SignatureMenu(object):
    def __init__(self, parser):
        required_group = parser.add_argument_group('Required arguments')
        required_group.add_argument("--infile", "-i", help="Input file in VCF or MAF format", type=argparse.FileType('r'))
        required_group.add_argument('--genome', "-g", help="Location of genome assembly file in 2bit format", type=str)
        required_group.add_argument("--signatures", "-s", choices=[5, 10, 30, 49, 53], help="Collection of signatures to use", type=int)

        optional_group = parser.add_argument_group('Optional arguments')
        optional_group.add_argument('--input-format', "-f", help="Input format: MAF, VCF", type=str, choices=['MAF', 'VCF', 'TCGI'], default='MAF')
        optional_group.add_argument(
            '--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout,
            help="Name of output file, will be generated in TSV format")

        # for backwards compatibility with 0.8.X add a hidden action that would just take anything as a valid input
        optional_group.add_argument('action', nargs="?", metavar="")

        advanced_group = parser.add_argument_group('Advanced arguments')
        advanced_group.add_argument('--method', "-m", help="Method defines the function minimized in the optimization procedure", type=str, default='MLEZ', nargs='?')
        advanced_group.add_argument('--no-unexplained-variance', "-U", help="Do not account for unexplained variance (non-context dependent mutational processes and unknown signatures)", action='store_false')
        advanced_group.add_argument('--bootstrap', "-b", help="Use the bootstrap to calculate confidence intervals", action='store_true')
        self.parser = parser

    def identify(self, args):
        if not args.infile:
            logger.warning("Provide input file in VCF or MAF format (-i) and a corresponding genome assembly (-g)")
            return
        if not args.genome:
            logger.warning(genome_error_message)
            return
        if not args.signatures:
            logger.warning("Set of signatures required. Use 5 and 10 for MUTAGENE-5 and MUTAGENE-10. Use 30 for COSMIC-30")
            return

        if args.method.lower() not in IDENTIFY_MIN_FUNCTIONS:
            logger.warning('Unknown method provided')
            return

        method = args.method

        # mutations, processing_stats = read_VCF_profile(args.infile, asm=args.genome)
        # mutations, processing_stats = read_auto_profile(args.infile, fmt=args.input_format, asm=args.genome)
        W, signature_names = read_signatures(int(args.signatures))

        if args.input_format == 'TCGI':
            mutations, mutations_with_context, processing_stats = read_TCGI_with_context_window(args.infile, args.genome, window_size=1)
            samples_profiles = get_multisample_mutational_profile(mutations, counts=True)

            samples_results = {}
            for sample, profile in samples_profiles.items():
                _, _, results = decompose_mutational_profile_counts(
                    profile,
                    (W, signature_names),
                    method,
                    others_threshold=0.0)
                samples_results[sample] = results
            write_multisample_decomposition(args.outfile, samples_results, signature_names)
        if args.input_format == 'MAF':
            mutations, mutations_with_context, processing_stats = read_MAF_with_context_window(args.infile, args.genome, window_size=1)
            samples_profiles = get_multisample_mutational_profile(mutations, counts=True)

            samples_results = {}
            for sample, profile in samples_profiles.items():
                _, _, results = decompose_mutational_profile_counts(
                    profile,
                    (W, signature_names),
                    method,
                    others_threshold=0.0)
                samples_results[sample] = results
            write_multisample_decomposition(args.outfile, samples_results, signature_names)
        elif args.input_format == 'VCF':
            mutations, processing_stats = read_auto_profile(args.infile, fmt=args.input_format, asm=args.genome)
            profile = get_mutational_profile(mutations, counts=True)
            if not args.bootstrap:
                _, _, results = decompose_mutational_profile_counts(
                    profile,
                    (W, signature_names),
                    method,
                    others_threshold=0.0,
                    enable_dummy=args.no_unexplained_variance)
                write_decomposition(args.outfile, results, signature_names, 'VCF')
            else:
                bootstrap_results = []
                for resampled_profile in generate_resampled_profiles(profile, 100):
                    _, _, results = decompose_mutational_profile_counts(
                        resampled_profile,
                        (W, signature_names),
                        method,
                        others_threshold=0.0,
                        enable_dummy=args.no_unexplained_variance)
                    bootstrap_results.append(results)
                write_bootstrap_decomposition(args.outfile, bootstrap_results, signature_names, 'VCF')

    def callback(self, args):
        self.identify(args)
