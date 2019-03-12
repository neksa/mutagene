import argparse
import signal
import sys

from pathlib import Path

import mutagene

# from mutagene.io import format_profile
from mutagene.io import read_profile_file, read_signatures
from mutagene.io import read_VCF_profile, read_MAF_with_genomic_context, get_mutational_profile, write_decomposition
from mutagene.io import read_cohort_mutations_from_tar
from mutagene.io import fetch_genome, fetch_cohorts
from mutagene.io import read_cohort_size_from_profile_file, list_cohorts_in_tar
from mutagene.profile import calc_profile
from mutagene.mutability import rank
from mutagene.version import __version__

import logging
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        prog="MutaGene",
        description='MutaGene version {} - Analysis of mutational processes and driver mutations'.format(mutagene.__version__),
        # usage="%(prog)s [options]",
        # formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('cmd', choices=['fetch_cohorts', 'fetch_genome', 'calc_profile', 'rank'])
    # parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                     help='an integer for the accumulator')
    # parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                     const=sum, default=max,
    #                     help='sum the integers (default: find the max)')
    parser.add_argument("--infile", "-i", nargs='*', help="Input file format", type=argparse.FileType('r'))
    parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    # parser.add_argument('--signatures', type=str)
    # parser.add_argument('--motifs', type=str)
    parser.add_argument('--genome', "-g", help="Location of genome assembly file", type=str)

    # parser.add_argument('--mode', "-m", help="n (nucleotide) or aa (amino acid) mode", type=str, default="aa")

    parser.add_argument('--cohorts-file', type=str, help="Location of tar.gz container or directory for cohorts", default="cohorts.tar.gz")
    parser.add_argument('--cohort', "-c", type=str, help="Name of cohort with observed mutations")
    parser.add_argument('--nsamples', "-n", type=int, help="Cohort size, overrides the value in profile")
    parser.add_argument('--profile', "-p", help="profile to calculate mutability, may also describe cohort size", type=str)

    # parser.add_argument('--threshold', "-t", help="B-score thresholds for drivers", type=float)
    # parser.add_argument('--mut-rate', "-r", help="Mutation rate overrides the rate inferred from profile", type=float)
    parser.add_argument('-v', '--verbose', action='count', default=0)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    args = parser.parse_args()

    # LOGGER INIT
    levels = [logging.WARNING, logging.INFO, logging.DEBUG]
    level = levels[min(len(levels) - 1, args.verbose)]  # capped to number of levels

    logging.basicConfig(level=level,
                        format="%(levelname)s %(message)s")

    genome_error_message = 'requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more'

    if args.cmd == 'fetch_genome':
        if not args.genome:
            logger.warning(genome_error_message)
            return
        fetch_genome(args.genome)
        logger.info("Twobit file saved to current directory")

    if args.cmd == 'fetch_cohorts':
        fetch_cohorts()
        logger.info("cohorts.tar.gz saved to current directory")

    if args.cmd == 'calc_profile':
        if not args.infile:
            logger.warning("Provide input file in VCF or MAF format (-i) and a corresponding genome assembly (-g)")
            return
        if not args.genome:
            logger.warning(genome_error_message)
            return
        calc_profile(args.infile, args.outfile, args.genome)

    if args.cmd == 'rank':
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


def identify_motifs(input_file, output_file, signatures, genome):
    from mutagene.motifs import decompose_mutational_profile_counts

    mutations, processing_stats = read_VCF_profile(input_file, asm=genome)
    W, signature_names = read_signatures(signatures)
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


def signal_handler(signal, frame):
    # logger.warning('Interrupted')
    sys.exit(0)


if __name__ == '__main__':
    signal.signal(signal.SIGINT, signal_handler)
    main()
