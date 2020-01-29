import sys
import os
# from pathlib import Path

# from mutagene.mutability.mutability import rank, THRESHOLD_DRIVER, THRESHOLD_PASSENGER
# from mutagene.io.cohorts import read_cohort_mutations_from_tar
# from mutagene.io.cohorts import read_cohort_size_from_profile_file, list_cohorts_in_tar
# from mutagene.io.profile import read_profile_file
# from mutagene.io.protein_mutations_MAF import read_MAF_with_genomic_context

from mutagene.benchmark.generate_benchmark import (
    gen_benchmark_2combinations,
    run_benchmark_2combinations,
    run_benchmark_2combinations_deconstruct_sigs,
    aggregate_benchmarks
)

from mutagene.benchmark.multiple_benchmark import (
    multiple_benchmark,
    multiple_benchmark_run,
    aggregate_multiple_benchmarks
)

from mutagene.io.profile import read_signatures

import logging
logger = logging.getLogger(__name__)


class BenchmarkMenu(object):
    def __init__(self, parser):
        required_group = parser.add_argument_group('Required arguments')
        required_group.add_argument("--mode", "-m", choices=[
            'pairwise_gen',
            'pairwise_run',
            'pairwise_run_ds',
            'multiple_gen',
            'multiple_run',
            'multiple_run_ds',
            'aggregate'], help="Multiple or pairwise mode, etc", type=str)
        required_group.add_argument("--signatures", "-i", nargs='*', help="Signatures (5, 10, 30, ...), default 30", type=str, default=["30"])

        # dirname = os.path.dirname(os.path.realpath(__file__))
        # default_root = dirname + "/../data/benchmark"
        default_root = "data/benchmark"
        default_root = os.path.normpath(default_root)

        required_group.add_argument(
            "--root",
            help="path to benchmark data, default {}".format(default_root),
            type=str,
            default=default_root)

        # required_group.add_argument('--genome', "-g", help="Location of genome assembly file in 2bit format", type=str)
        self.parser = parser
        pass

    def callback(self, args):
        if args.mode.startswith("pair"):
            for i in args.signatures:
                i = int(i)
                W, signature_names = read_signatures(i)

                if args.mode.endswith('gen'):
                    gen_benchmark_2combinations(args.root, signature_names, W)
                elif args.mode.endswith('run'):
                    run_benchmark_2combinations(args.root, i, signature_names, W, force=True)
                elif args.mode.endswith('run_ds'):
                    run_benchmark_2combinations_deconstruct_sigs(args.root, i, signature_names, W, force=True)

        elif args.mode.startswith('multiple'):
            if args.mode.endswith('gen'):
                pass
            elif args.mode.endswith('run'):
                multiple_benchmark_run(i, signature_names, W, force=True)
            elif args.mode.endswith('run_ds'):
                pass
            multiple_benchmark()
            aggregate_multiple_benchmarks()

        elif args.mode == 'aggregate':
            aggregate_benchmarks(args.root)

        else:
            print("Unknown benchmark action mode")
            self.parser.print_usage()
            sys.exit(1)
