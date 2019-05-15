import argparse
import signal
import sys

from pathlib import Path
import mutagene
from mutagene.version import __version__
from mutagene.io.fetch import fetch_genome, fetch_cohorts, fetch_examples

from mutagene.io.profile import format_profile
from mutagene.io.profile import read_profile_file, read_signatures
from mutagene.io.mutations_profile import read_VCF_profile
from mutagene.io.protein_mutations_MAF import read_MAF_with_genomic_context
from mutagene.profiles.profile import get_mutational_profile
# from mutagene.io.decomposition import write_decomposition
from mutagene.io.cohorts import read_cohort_mutations_from_tar
from mutagene.io.cohorts import read_cohort_size_from_profile_file, list_cohorts_in_tar
from mutagene.mutability.mutability import rank

import logging
logger = logging.getLogger(__name__)


genome_error_message = 'requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more'


class FetchMenu(object):
    def __init__(self, parser):
        parser.description = 'Download data from remote repositories and API'
        # parser.add_argument('resource', choices=['cohorts', 'genome', 'GDC', 'MSKC', 'ICGC'])
        subparsers = parser.add_subparsers(
            dest='resource',
            title='subcommands',
            description='Choose data source',
            help='additional help available for subcommands')
        examples_parser = subparsers.add_parser('examples')
        cohorts_parser = subparsers.add_parser('cohorts')
        cohorts_parser.add_argument('--cohort', choices=('COSMIC', 'GDC', 'MSKC', 'ICGC'))
        genome_parser = subparsers.add_parser('genome')
        genome_parser.add_argument('--genome', '-g', type=str, help='hg38, hg19, mm10 according to UCSC genome browser nomenclature')

    @classmethod
    def examples(cls, args):
        # print('Cohorts', args)
        if args.resource == 'examples':
            fetch_examples()
            logger.info("Example files saved to current directory")

    @classmethod
    def cohorts(cls, args):
        # print('Cohorts', args)
        if args.resource == 'cohorts':
            fetch_cohorts()
            logger.info("cohorts.tar.gz saved to current directory")

    @classmethod
    def genome(cls, args):
        # print('Genome', args)
        if args.resource == 'genome':
            if not args.genome:
                logger.warning(genome_error_message)
                return
            fetch_genome(args.genome)
            logger.info("Twobit file saved to current directory")

    @classmethod
    def callback(cls, args):
        if not args.resource:
            logger.warning('No resource specified')
            sys.exit(1)

        print('FetchMenu', args.resource)
        getattr(cls, args.resource)(args)

################### PROFILE ####################
from mutagene.profiles.profile import calc_profile
class ProfileMenu(object):
    def __init__(self, parser):
        parser.add_argument('action', choices=['calculate'])  # , 'plot', 'fingerprint', 'compare'])
        parser.add_argument("--infile", "-i", nargs='*', help="Input file format", type=argparse.FileType('r'))
        parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
        parser.add_argument('--genome', "-g", help="Location of genome assembly file", type=str)

    @classmethod
    def callback(cls, args):
        print('ProfileMenu', args.action)
        getattr(cls, args.action)(args)

    @classmethod
    def calculate(cls, args):
        print("Calculating...")

        if not args.infile:
            logger.warning("Provide input file in VCF or MAF format (-i) and a corresponding genome assembly (-g)")
            return
        if not args.genome:
            logger.warning(genome_error_message)
            return
        calc_profile(args.infile, args.outfile, args.genome)

################################################
################################################


class MotifMenu(object):
    def __init__(self, parser):
        parser.add_argument('action', choices=['search'])
        parser.add_argument("--motif", "-m", help="Motif to search for, use the 'R[C>T]GY' syntax for the motif. Use quotes", type=str)
        parser.add_argument("--infile", "-i", help="Input file in MAF or VCF format with one or multiple samples", type=argparse.FileType('r'))
        parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
        parser.add_argument('--genome', "-g", help="Location of genome assembly file in 2bit format", type=str)

    @classmethod
    def search(cls, args):
        from mutagene.io.context_window import read_MAF_with_context_window
        from mutagene.motifs.motifs import identify_motifs
        from mutagene.io.motifs import write_motif_matches

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

        mutations, mutations_with_context, processing_stats = read_MAF_with_context_window(args.infile, args.genome)
        matching_motifs = identify_motifs(mutations_with_context, custom_motif) if mutations_with_context is not None else []
        if len(matching_motifs) == 0:
            logger.warning("No significant motif matches found")
        else:
            write_motif_matches(args.outfile, matching_motifs)

    @classmethod
    def callback(cls, args):
        # print('MotifMenu', args.action)
        getattr(cls, args.action)(args)


class SignatureMenu(object):
    def __init__(self, parser):
        parser.add_argument('action', choices=['identify', ])  # 'new'])
        parser.add_argument("--signatures", "-s", choices=[5, 10, 30], help="Collection of signatures to use", type=int)
        parser.add_argument("--infile", "-i", help="Input file in VCF or MAF format", type=argparse.FileType('r'))
        parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
        parser.add_argument('--genome', "-g", help="Location of genome assembly file in 2bit format", type=str)

    @classmethod
    def identify(cls, args):
        # def identify_signatures(input_file, output_file, signatures, genome):
        from mutagene.io.mutations_profile import read_auto_profile
        from mutagene.signatures.identify import decompose_mutational_profile_counts
        from mutagene.io.decomposition import write_decomposition

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


class RankMenu(object):
    def __init__(self, parser):
        parser.add_argument("--infile", "-i", nargs='*', help="Input file in MAF format", type=argparse.FileType('r'))
        parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
        parser.add_argument('--genome', "-g", help="Location of genome assembly file", type=str)

        # parser.add_argument('--mode', "-m", help="n (nucleotide) or aa (amino acid) mode", type=str, default="aa")

        parser.add_argument('--cohorts-file', type=str, help="Location of tar.gz container or directory for cohorts", default="cohorts.tar.gz")
        parser.add_argument('--cohort', "-c", type=str, help="Name of cohort with observed mutations")
        parser.add_argument('--nsamples', "-n", type=int, help="Cohort size, overrides the value in profile")
        parser.add_argument('--profile', "-p", help="profile to calculate mutability, may also describe cohort size", type=str)

    @classmethod
    def callback(cls, args):
        # print('RankMenu', args.infile)
        # if not args.resource:
        #     logger.warning('No resource specified')
        #     sys.exit(1)

        # getattr(cls, args.resource)(args)

        # if args.mode == 'na':
        #     logger.warning('DNA mutations ranking not supported')
        #     return

        if not args.genome:
            logger.warning(genome_error_message)
            return

        if not Path(args.cohorts_file).is_file():
            logger.warning("Cohorts file missing. Download with \"mutagene fetch_cohorts\"")
            return

        if args.cohort and Path(args.cohorts_file).is_file():
            profile, cohort_size, cohort_aa_mutations, cohort_na_mutations = read_cohort_mutations_from_tar(args.cohorts_file, args.cohort)
            logger.info('Profile and cohort size loaded from precalculated cohorts N=' + str(cohort_size))
        else:
            logger.warning('Cohort required')

            if args.cohorts_file:
                logger.warning("List of available cohorts:\n" + list_cohorts_in_tar(args.cohorts_file))
            return

        if args.profile:
            profile = read_profile_file(args.profile)
            if profile:
                logger.info('Profile overridden')
            else:
                return
            cohort_size_new = read_cohort_size_from_profile_file(args.profile)
            if cohort_size_new:
                cohort_size = cohort_size_new
                logger.info('Cohort size loaded from profile N=' + str(cohort_size))

        if args.nsamples:
            cohort_size = args.nsamples
            logger.info('Cohort size overridden N=' + str(cohort_size))

        if isinstance(args.infile, list):
            if len(args.infile) > 1:
                logger.info('Multiple input files provided')
            mutations_to_rank = []
            processing_stats = {'loaded': 0, 'skipped': 0}
            for infile in args.infile:
                mut, stats = read_MAF_with_genomic_context(infile, args.genome)
                mutations_to_rank.extend(mut)
                processing_stats['loaded'] += stats['loaded']
                processing_stats['skipped'] += stats['skipped']
        else:
            mutations_to_rank, processing_stats = read_MAF_with_genomic_context(args.infile, args.genome)

        if not len(mutations_to_rank):
            logger.warning('No mutations to rank')
            return

        msg = "Loaded {} mutations".format(processing_stats['loaded'])
        if processing_stats['skipped'] > 0:
            msg += " skipped {} mutations".format(processing_stats['skipped'])
        logger.info(msg)

        rank(mutations_to_rank, args.outfile, profile, cohort_aa_mutations, cohort_size)


"""

# parser.add_argument('integers', metavar='N', type=int, nargs='+',
#                     help='an integer for the accumulator')
# parser.add_argument('--sum', dest='accumulate', action='store_const',
#                     const=sum, default=max,
#                     help='sum the integers (default: find the max)')
parser.add_argument("--infile", "-i", nargs='*', help="Input file format", type=argparse.FileType('r'))
parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
# parser.add_argument('--signatures', type=str)
# parser.add_argument('--motifs', type=str)

# parser.add_argument('--threshold', "-t", help="B-score thresholds for drivers", type=float)
# parser.add_argument('--mut-rate', "-r", help="Mutation rate overrides the rate inferred from profile", type=float)

"""


class MutaGeneApp(object):
    def __init__(self):
        signal.signal(signal.SIGINT, self.signal_handler)

        parser = argparse.ArgumentParser(
            prog="MutaGene",
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



    # if args.cmd == 'rank':
    #     # if args.mode == 'na':
    #     #     logger.warning('DNA mutations ranking not supported')
    #     #     return

    #     if not args.genome:
    #         logger.warning(genome_error_message)
    #         return

    #     if not Path(args.cohorts_file).is_file():
    #         logger.warning("Cohorts file missing. Download with \"mutagene fetch_cohorts\"")
    #         return

    #     if args.cohort and Path(args.cohorts_file).is_file():
    #         profile, cohort_size, cohort_aa_mutations, cohort_na_mutations = read_cohort_mutations_from_tar(args.cohorts_file, args.cohort)
    #         logger.info('Profile and cohort size loaded from precalculated cohorts N=' + str(cohort_size))
    #     else:
    #         logger.warning('Cohort required')

    #         if args.cohorts_file:
    #             logger.warning("List of available cohorts:\n" + list_cohorts_in_tar(args.cohorts_file))
    #         return

    #     if args.profile:
    #         profile = read_profile_file(args.profile)
    #         if profile:
    #             logger.info('Profile overridden')
    #         else:
    #             return
    #         cohort_size_new = read_cohort_size_from_profile_file(args.profile)
    #         if cohort_size_new:
    #             cohort_size = cohort_size_new
    #             logger.info('Cohort size loaded from profile N=' + str(cohort_size))

    #     if args.nsamples:
    #         cohort_size = args.nsamples
    #         logger.info('Cohort size overridden N=' + str(cohort_size))

    #     if isinstance(args.infile, list):
    #         logger.info('Multiple input files provided')
    #         mutations_to_rank = []
    #         processing_stats = {'loaded': 0, 'skipped': 0}
    #         for infile in args.infile:
    #             mut, stats = read_MAF_with_genomic_context(infile, args.genome)
    #             mutations_to_rank.extend(mut)
    #             processing_stats['loaded'] += stats['loaded']
    #             processing_stats['skipped'] += stats['skipped']
    #     else:
    #         mutations_to_rank, processing_stats = read_MAF_with_genomic_context(args.infile, args.genome)

    #     if not len(mutations_to_rank):
    #         logger.warning('No mutations to rank')
    #         return

    #     msg = "Loaded {} mutations".format(processing_stats['loaded'])
    #     if processing_stats['skipped'] > 0:
    #         msg += " skipped {} mutations".format(processing_stats['skipped'])
    #     logger.info(msg)

    #     rank(mutations_to_rank, args.outfile, profile, cohort_aa_mutations, cohort_size)

    """
    if args.cmd == 'benchmark':
        benchmark()

    if args.cmd == 'identify':
        if not args.genome:
            logger.warning(genome_error_message)
            return
        if args.signatures:
            identify_signatures(args.infile, args.outfile, args.signatures, args.genome)
        elif args.motifs:
            identify_motifs(args.infile, args.outfile, args.signatures, args.genome)
    """


"""
def identify_signatures(input_file, output_file, signatures, genome):
    from mutagene.identify import decompose_mutational_profile_counts

    mutations, processing_stats = read_VCF_profile(input_file, asm=genome)
    W, signature_names = read_signatures(int(signatures))
    # print(W, signature_names)
    profile = get_mutational_profile(mutations, counts=True)

    _, _, results = decompose_mutational_profile_counts(
        profile,
        (W, signature_names),
        'MLEZ',
        debug=False,
        others_threshold=0.0)
    write_decomposition(output_file, results, signature_names)

# def gdc():
#     from .gdc import gdc_read_file
#     print(gdc_read_file())


def global_optimization(input_file):
    from .identify import decompose_mutational_profile_counts

    profile = read_profile_file(input_file)
    W, signature_names = read_signatures(30)
    _, _, results = decompose_mutational_profile_counts(
        profile,
        (W, signature_names),
        'MLEZ-GLOB',
        debug=False,
        others_threshold=0.0)
    print(results)


def benchmark():
    from .generate_benchmark import gen_benchmark_2combinations
    from .generate_benchmark import run_benchmark_2combinations, run_benchmark_2combinations_deconstruct_sigs
    from .generate_benchmark import aggregate_benchmarks
    from .multiple_benchmark import multiple_benchmark, multiple_benchmark_run, aggregate_multiple_benchmarks

    ##### multiple
    # multiple_benchmark()

    ##### pairwise
    # for i in [5, 10, 30]:
    for i in [30, ]:
        W, signature_names = read_signatures(i)
        # gen_benchmark_2combinations(signature_names, W)
        # run_benchmark_2combinations_deconstruct_sigs(i, signature_names, W, force=True)
        # multiple_benchmark_run(i, signature_names, W, force=True)
        # run_benchmark_2combinations(i, signature_names, W, force=True)

    aggregate_multiple_benchmarks()
    # aggregate_benchmarks()
"""


if __name__ == '__main__':
    MutaGeneApp()
