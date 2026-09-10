import logging
import os
import shutil
import sys

# from pathlib import Path
# from mutagene.mutability.mutability import rank, THRESHOLD_DRIVER, THRESHOLD_PASSENGER
# from mutagene.io.cohorts import read_cohort_mutations_from_tar
# from mutagene.io.cohorts import read_cohort_size_from_profile_file, list_cohorts_in_tar
# from mutagene.io.profile import read_profile_file
# from mutagene.io.protein_mutations_MAF import read_MAF_with_genomic_context
from mutagene.benchmark.generate_benchmark import (
    aggregate_benchmarks,
    gen_benchmark_2combinations,
    run_benchmark_2combinations,
    run_benchmark_2combinations_deconstruct_sigs,
)
from mutagene.benchmark.multiple_benchmark import (
    aggregate_multiple_benchmarks,
    multiple_benchmark,
    multiple_benchmark_run,
)
from mutagene.io.profile import read_signatures

logger = logging.getLogger(__name__)


class BenchmarkMenu:
    def __init__(self, parser):
        required_group = parser.add_argument_group("Required arguments")
        required_group.add_argument(
            "--mode",
            "-m",
            choices=[
                "pairwise_gen",
                "pairwise_run",
                "pairwise_run_ds",
                "multiple_gen",
                "multiple_run",
                "multiple_run_ds",
                "aggregate",
            ],
            help="Multiple or pairwise mode, etc",
            type=str,
        )
        required_group.add_argument(
            "--signatures",
            "-i",
            nargs="*",
            help="Signatures (5, 10, 30, ...), default 30",
            type=str,
            default=["30"],
        )

        # dirname = os.path.dirname(os.path.realpath(__file__))
        # default_root = dirname + "/../data/benchmark"
        default_root = "data/benchmark"
        default_root = os.path.normpath(default_root)

        required_group.add_argument(
            "--root",
            help=f"path to benchmark data, default {default_root}",
            type=str,
            default=default_root,
        )

        optional = parser.add_argument_group("Optional arguments")
        optional.add_argument(
            "--replicates",
            type=int,
            default=100,
            help="Number of synthetic samples to generate (multiple_gen), default 100",
        )
        optional.add_argument(
            "--processes",
            type=int,
            default=10,
            help="Worker processes, default 10. Use 1 to see errors directly",
        )

        # required_group.add_argument('--genome', "-g", help="Location of genome assembly file in 2bit format", type=str)
        self.parser = parser
        pass

    def callback(self, args):
        # read_signatures keys on strings; int() made every pairwise mode fail
        # with "Unknown signature set: 5" against a list that visibly contains 5.
        signature_sets = [str(s) for s in args.signatures]

        if args.mode.endswith("_ds") and not shutil.which("Rscript"):
            logger.error(
                "The deconstructSigs modes need R with the deconstructSigs package, "
                "and Rscript is not on PATH"
            )
            sys.exit(1)

        if args.mode.startswith("pairwise"):
            for name in signature_sets:
                W, signature_names = read_signatures(name)

                if args.mode == "pairwise_gen":
                    gen_benchmark_2combinations(args.root, signature_names, W)
                elif args.mode == "pairwise_run":
                    run_benchmark_2combinations(args.root, name, signature_names, W, force=True)
                elif args.mode == "pairwise_run_ds":
                    run_benchmark_2combinations_deconstruct_sigs(
                        args.root, name, signature_names, W, force=True
                    )

        elif args.mode.startswith("multiple"):
            # These used i, signature_names and W, which are only ever assigned
            # in the pairwise branch, so every multiple mode raised NameError.
            if args.mode == "multiple_gen":
                multiple_benchmark(
                    data_root=args.root,
                    signature_sets=signature_sets,
                    replicates=args.replicates,
                    processes=args.processes,
                )
            elif args.mode in ("multiple_run", "multiple_run_ds"):
                for name in signature_sets:
                    W, signature_names = read_signatures(name)
                    multiple_benchmark_run(
                        name,
                        signature_names,
                        W,
                        force=True,
                        data_root=args.root,
                        processes=args.processes,
                    )
                aggregate_multiple_benchmarks(data_root=args.root)

        elif args.mode == "aggregate":
            aggregate_benchmarks(args.root)

        else:
            logger.error(f"Unknown benchmark mode: {args.mode}")
            self.parser.print_usage()
            sys.exit(1)
