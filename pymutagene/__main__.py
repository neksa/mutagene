import argparse
import sys

from .io import read_profile, format_profile, read_signatures
from .io import read_VCF, get_mutational_profile, write_decomposition


def main():
    parser = argparse.ArgumentParser(description='Analysis of mutational processes with PyMutaGene package')
    # parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                     help='an integer for the accumulator')
    # parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                     const=sum, default=max,
    #                     help='sum the integers (default: find the max)')
    parser.add_argument('cmd', choices=['benchmark', 'identify', 'gdc', 'global'])
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'))
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()

    if args.cmd == 'benchmark':
        benchmark()

    if args.cmd == 'identify':
        identify(args.infile, args.outfile)

    if args.cmd == 'gdc':
        gdc()

    if args.cmd == 'global':
        global_optimization(args.infile)


def gdc():
    from .gdc import gdc_read_file
    print(gdc_read_file())


def global_optimization(input_file):
    from .identify import decompose_mutational_profile_counts

    profile = read_profile(input_file)
    W, signature_names = read_signatures(30)
    _, _, results = decompose_mutational_profile_counts(
        profile,
        (W, signature_names),
        'MLEZ-GLOB',
        debug=False,
        others_threshold=0.0)
    print(results)


def identify(input_file, output_file, signatures=30):
    from .identify import decompose_mutational_profile_counts

    mutations, processing_stats = read_VCF(input_file)
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
