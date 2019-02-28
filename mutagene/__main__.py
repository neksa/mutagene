import logging
import argparse
import sys

import mutagene

from mutagene.io import read_profile_file, format_profile, read_signatures
from mutagene.io import read_VCF_profile, read_MAF_with_genomic_context, get_mutational_profile, write_decomposition
from mutagene.io import read_cohort_mutations_from_tar

from mutagene.profile import calc_profile

from mutagene.mutability import rank


def main():
    logging.basicConfig(level=logging.DEBUG)

    parser = argparse.ArgumentParser(
        description='MutaGene version {} - Analysis of mutational processes and driver mutations'.format(mutagene.__version__),
        # usage="%(prog)s [options]",
        # formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('cmd', choices=['identify', 'rank', 'calc_freq', 'calc_profile', 'benchmark', ])
    # parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                     help='an integer for the accumulator')
    # parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                     const=sum, default=max,
    #                     help='sum the integers (default: find the max)')
    parser.add_argument("--infile", "-i", nargs='?', help="Input file format", type=argparse.FileType('r'))
    parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--signatures', type=str)
    parser.add_argument('--motifs', type=str)
    # parser.add_argument('--assembly', "-a", type=str)
    parser.add_argument('--cohort', "-c", type=str, help="Name of cohort with observed mutations")
    parser.add_argument('--nsamples', "-n", type=str, help="Cohort size, overrides the value in profile")
    parser.add_argument('--cohorts-file', type=str, help="Location of tar.gz container or directory for cohorts", default="cohorts.tar.gz")
    parser.add_argument('--genome', "-g", help="Location of genome assembly file", type=str)
    # parser.add_argument('--normalization', "-n", help="Nucleotide composition normalization: exome_xy, genome_xy(default)", type=str)
    parser.add_argument('--mode', "-m", help="n (nucleotide) or aa (amino acid) mode", type=str, default="aa")
    parser.add_argument('--profile', "-p", help="profile to calculate mutability, may also describe cohort size", type=str)
    parser.add_argument('--threshold', "-t", help="B-score thresholds for drivers", type=float)
    parser.add_argument('--mut-rate', "-r", help="Mutation rate overrides the rate inferred from profile", type=float)
    args = parser.parse_args()

    if args.cmd == 'benchmark':
        benchmark()

    if args.cmd == 'calc_freq':
        calc_freq()

    if args.cmd == 'calc_profile':
        calc_profile(args.infile, args.outfile, args.genome)

    if args.cmd == 'identify':
        if args.signatures:
            identify_signatures(args.infile, args.outfile, args.signatures, args.genome)
        elif args.motifs:
            identify_motifs(args.infile, args.outfile, args.signatures, args.genome)

    if args.cmd == 'rank':
        if args.mode == 'na':
            logging.warning('DNA mutations ranking not supported')
            return

        if args.cohort and args.cohorts_file:
            profile, cohort_size, cohort_aa_mutations, cohort_na_mutations = read_cohort_mutations_from_tar(args.cohorts_file, args.cohort)
        else:
            logging.warning('Cohort required')
            return

        mutations_to_rank = read_MAF_with_genomic_context(args.infile, args.genome)

        results = rank(mutations_to_rank, args.outfile, profile, cohort_aa_mutations, cohort_size)
        import pprint
        pprint.pprint(results)
        
    # if args.cmd == 'gdc':
    #     gdc()

    # if args.cmd == 'global':
    #     global_optimization(args.infile)


def gdc():
    from .gdc import gdc_read_file
    print(gdc_read_file())


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


if __name__ == '__main__':
    main()
