import argparse

def main():
    parser = argparse.ArgumentParser(description='Analysis of mutational processes with PyMutaGene package')
    # parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                     help='an integer for the accumulator')
    # parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                     const=sum, default=max,
    #                     help='sum the integers (default: find the max)')
    parser.add_argument('cmd', choices=['benchmark', 'run'])
    args = parser.parse_args()
    if args.cmd == 'benchmark':
        benchmark()


def benchmark():
    import numpy as np
    from .benchmark import benchmark_simulated, benchmark_2combinations
    from .io import read_profile, format_profile

    W = []
    signature_names = []
    for i in range(5):
        fname = "data/signatures/A_{}.profile".format(i + 1)
        profile = read_profile(fname)
        W.append(profile)
        signature_names.append("{}".format(i + 1))

    W = np.array(W)
    # print(W.shape)
    # print(signature_names)
    # benchmark_simulated("/tmp/test_benchmark.txt", signature_names, W)

    benchmark_2combinations("/tmp/test_benchmark2.txt", signature_names, W)


if __name__ == '__main__':
    main()
