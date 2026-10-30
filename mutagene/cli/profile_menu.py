import argparse
import logging
import sys

from mutagene.cli.genome_resolver import GenomeNotFoundError, resolve_genome
from mutagene.io.text import open_text_input
from mutagene.io.variant_filter import add_filter_argument
from mutagene.profiles.profile import calc_profile

logger = logging.getLogger(__name__)
genome_error_message = """requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.soe.ucsc.edu/downloads.html for more
                          Use mutagene fetch to download genome assemblies"""


class ProfileMenu:
    def __init__(self, parser):
        parser.add_argument(
            "--infile", "-i", nargs="*", help="Input file format", type=open_text_input
        )
        parser.add_argument(
            "--outfile",
            "-o",
            nargs="?",
            type=argparse.FileType("w"),
            default=sys.stdout,
            help="Name of output file, will be generated in TSV format",
        )
        parser.add_argument("--genome", "-g", help="Location of genome assembly file", type=str)
        parser.add_argument(
            "--input-format", "-f", help="Input format: auto, MAF, VCF", type=str, default="auto"
        )
        parser.add_argument(
            "--per-sample",
            action="store_true",
            help="Write a labelled SBS96 count matrix with one column per sample",
        )

        # for backwards compatibility with 0.8.X add a hidden action that would just take anything as a valid input
        add_filter_argument(parser)

        parser.add_argument("action", nargs="?", metavar="")

    def callback(self, args):
        # print('ProfileMenu', args.action)
        self.calculate(args)

    def calculate(self, args):
        # print("Calculating...")
        if not args.infile:
            logger.warning(
                "Provide input file in VCF or MAF format (-i) and a corresponding genome assembly (-g)"
            )
            sys.exit(1)
        if not args.genome:
            logger.warning(genome_error_message)
            sys.exit(1)
        try:
            args.genome = resolve_genome(args.genome)
        except GenomeNotFoundError as e:
            logger.error(str(e))
            sys.exit(1)

        try:
            if args.per_sample:
                from mutagene.io.context_window import read_mutations
                from mutagene.io.matrix import write_count_matrix
                from mutagene.io.text import detect_mutation_format
                from mutagene.profiles.profile import get_multisample_mutational_profile

                profiles = {}
                for stream in args.infile:
                    fmt = args.input_format.upper()
                    if fmt == "AUTO":
                        fmt = detect_mutation_format(stream)
                    mutations, _, _ = read_mutations(
                        fmt,
                        stream,
                        args.genome,
                        window_size=1,
                        keep_filtered=args.keep_filtered,
                        filter_column=args.filter_column,
                    )
                    current = get_multisample_mutational_profile(mutations, counts=True)
                    if profiles.keys() & current.keys():
                        raise ValueError(
                            "Sample names repeat across input files; supply unique sample names"
                        )
                    profiles.update(current)
                if not profiles:
                    raise ValueError("No mutations could be read for the sample matrix")
                write_count_matrix(args.outfile, profiles)
                return
            stats = calc_profile(
                args.infile,
                args.outfile,
                args.genome,
                args.input_format,
                keep_filtered=args.keep_filtered,
                filter_column=args.filter_column,
            )
        except (ValueError, OSError, EOFError) as e:
            # An unreadable input is the user's problem to fix, not a defect to
            # report, so say what is wrong instead of printing a traceback.
            # rank and signature already did this; profile did not.
            logger.error(
                f"Could not read the input in {args.input_format} format: {e}\n"
                "Check the file, or name a different format with --input-format (-f)"
            )
            if logger.root.level == logging.DEBUG:
                raise
            sys.exit(1)

        # An empty result is indistinguishable from a successful one to anything downstream, so a run that loaded nothing has to say so in its exit status.
        if not stats or not stats.get("loaded"):
            logger.error(
                "No mutations could be read, so no profile was written. "
                "Check the input format and that the genome assembly matches it"
            )
            sys.exit(1)
