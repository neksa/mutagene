import argparse
import sys
import logging

from tqdm import tqdm

from mutagene.io.profile import read_signatures
from mutagene.profiles.profile import get_multisample_mutational_profile
from mutagene.profiles.profile import generate_resampled_profiles
from mutagene.io.context_window import read_mutations
from mutagene.signatures.identify import decompose_mutational_profile_counts
from mutagene.signatures.identify import IDENTIFY_MIN_FUNCTIONS
from mutagene.io.decomposition import write_decomposition


logger = logging.getLogger(__name__)

genome_error_message = """requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more
                          Use mutagene fetch to download genome assemblies"""


class SignatureMenu(object):
    def __init__(self, parser):
        required_group = parser.add_argument_group('Required arguments')
        required_group.add_argument("--infile", "-i", help="Input file in VCF or MAF format", type=argparse.FileType('r'))
        required_group.add_argument('--genome', "-g", help="Location of genome assembly file in 2bit format", type=str, default='hg19')
        required_group.add_argument(
            "--signatures", "-s",
            choices=["5", "10", "30", "49", "53", "MGA", "MGB", "COSMICv2", "COSMICv3", "KUCAB"],
            help="Collection of signatures to use", type=str, default="COSMICv2")

        optional_group = parser.add_argument_group('Optional arguments')
        optional_group.add_argument('--input-format', "-f", help="Input format: MAF, VCF", type=str, choices=['MAF', 'VCF', 'TCGI'], default='MAF')
        optional_group.add_argument(
            '--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout,
            help="Name of output file, will be generated in TSV format")

        # for backwards compatibility with 0.8.X add a hidden action that would just take anything as a valid input
        optional_group.add_argument('action', nargs="?", metavar="")

        advanced_group = parser.add_argument_group('Advanced arguments')
        advanced_group.add_argument('--method', "-m", help="Method defines the function minimized in the optimization procedure", type=str, default='MLE', nargs='?')
        advanced_group.add_argument('--no-unexplained-variance', "-U", help="Do not account for unexplained variance (non-context dependent mutational processes and unknown signatures)", action='store_false')
        advanced_group.add_argument('--mutations-threshold', "-t", help="Only report signatures with mutations above the threshold", type=int, default=0)
        advanced_group.add_argument('--keep-only', "-k", help="Keep only the signatures in the list, separated by commas e.g. 1,3,5", type=str, default=None)

        bootstrap_group = parser.add_argument_group('Bootstrap-specific arguments')
        bootstrap_group.add_argument('--bootstrap', "-b", help="Use the bootstrap to calculate confidence intervals", action='store_true')
        bootstrap_group.add_argument('--bootstrap-replicates', "-br", help="Number of bootstrap replicates", type=int, default=100)
        bootstrap_group.add_argument('--bootstrap-confidence-level', "-bcl", help="Confidence level", type=int, default=95)
        bootstrap_group.add_argument('--bootstrap-method', "-bm", help="Bootstrap method (t: t-distribution, p: percentile)", type=str, choices=['t', 'p'], default='p')
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

        if args.bootstrap_replicates < 10:
            logger.warning("Number of bootstrap replicates too low. Specify at least 10 replicates")
            return

        if args.bootstrap_confidence_level < 70:
            logger.warning("Specify confidence level of at least 70% and less than 99%")
            return

        only = None
        if args.keep_only is not None:
            only = args.keep_only.split(",")
            if len(only) < 1:
                logger.warning("List of signatures to keep for the analysis is empty")
                return
            logger.warning("We will only analyze signatures in this list: {}".format(", ".join(only)))

        if 'input_format' not in args:
            # guess format from file name
            name = args.infile.name.upper()
            if name.endswith("MAF"):
                args.input_format = "MAF"
            elif name.endswith("VCF"):
                args.input_format = "VCF"
            else:
                logger.warning("Input format was not specified. Assuming it is MAF")
                args.input_format = "MAF"

        W, signature_names = read_signatures(args.signatures, only=only)

        try:
            mutations, _, processing_stats = read_mutations(args.input_format, args.infile, args.genome, window_size=1)
        except Exception as e:
            e_message = getattr(e, 'message', repr(e))
            logger.warning(
                "Parsing {0} failed. "
                "Check that the input file is in {0} format "
                "or specify a different format using option -f \n"
                "{1}".format(args.input_format, e_message))

            if logger.root.level == logging.DEBUG:
                raise
            return

        samples_profiles = get_multisample_mutational_profile(mutations, counts=True)
        samples_results = {}

        for sample, profile in samples_profiles.items():
            _, _, results = decompose_mutational_profile_counts(
                profile,
                (W, signature_names),
                args.method,
                others_threshold=0.0,
                enable_dummy=args.no_unexplained_variance)
            samples_results[sample] = results

        if not args.bootstrap:
            write_decomposition(args.outfile, samples_results, signature_names, mutations_threshold=args.mutations_threshold)
        else:
            bootstrap_samples_results = {}
            for sample, profile in samples_profiles.items():
                bootstrap_results = []
                for resampled_profile in tqdm(generate_resampled_profiles(profile, args.bootstrap_replicates), total=args.bootstrap_replicates):
                    _, _, results = decompose_mutational_profile_counts(
                        resampled_profile,
                        (W, signature_names),
                        args.method,
                        others_threshold=0.0,
                        enable_dummy=args.no_unexplained_variance)
                    bootstrap_results.append(results)
                bootstrap_samples_results[sample] = bootstrap_results

            write_decomposition(
                args.outfile, samples_results, signature_names,
                mutations_threshold=args.mutations_threshold,
                bootstrap_method=args.bootstrap_method,
                profile=samples_profiles,
                bootstrap_results=bootstrap_samples_results,
                bootstrap_level=args.bootstrap_confidence_level)

    def callback(self, args):
        self.identify(args)
