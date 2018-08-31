import argparse


def main():
    parser = argparse.ArgumentParser(description='Analysis of mutational processes with PyMutaGene package')
    # parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                     help='an integer for the accumulator')
    # parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                     const=sum, default=max,
    #                     help='sum the integers (default: find the max)')
    parser.add_argument('cmd', choices=['benchmark'])
    args = parser.parse_args()

    if args.cmd == 'benchmark':
        benchmark()


def benchmark():
    from .io import read_profile, format_profile, read_signatures
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
