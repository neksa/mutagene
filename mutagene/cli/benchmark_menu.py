# import argparse
# import sys
import logging
# from pathlib import Path

# from mutagene.mutability.mutability import rank, THRESHOLD_DRIVER, THRESHOLD_PASSENGER
# from mutagene.io.cohorts import read_cohort_mutations_from_tar
# from mutagene.io.cohorts import read_cohort_size_from_profile_file, list_cohorts_in_tar
# from mutagene.io.profile import read_profile_file
# from mutagene.io.protein_mutations_MAF import read_MAF_with_genomic_context


logger = logging.getLogger(__name__)
genome_error_message = """requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more
                          Use mutagene fetch to download genome assemblies"""


class BenchmarkMenu(object):
    def __init__(self, parser):
        pass

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

    @classmethod
    def callback(cls, args):

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
