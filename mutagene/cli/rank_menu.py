import argparse
import sys
import logging
from pathlib import Path

from mutagene.mutability.mutability import rank, THRESHOLD_DRIVER, THRESHOLD_PASSENGER
from mutagene.io.cohorts import read_cohort_mutations_from_tar
from mutagene.io.cohorts import read_cohort_size_from_profile_file, list_cohorts_in_tar
from mutagene.io.profile import read_profile_file
from mutagene.io.context_window import read_mutations
from mutagene.profiles.profile import get_pooled_multisample_mutational_profile
from mutagene.io.protein_mutations_MAF import read_protein_mutations_MAF


logger = logging.getLogger(__name__)
genome_error_message = """requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more
                          Use mutagene fetch to download genome assemblies"""


class RankMenu(object):
    def __init__(self, parser):
        required_group = parser.add_argument_group('Required arguments')
        required_group.add_argument("--infile", "-i", help="Input file in MAF format", type=argparse.FileType('r'))
        required_group.add_argument('--genome', "-g", help="Location of genome assembly file in 2bit format", type=str, default='hg19')
        # parser.add_argument('--mode', "-m", help="n (nucleotide) or aa (amino acid) mode", type=str, default="aa")

        optional_group = parser.add_argument_group('Optional arguments')
        optional_group.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
        # Suppressed with issue #51
        optional_group.add_argument('--input-format', "-f", help=argparse.SUPPRESS, type=str, choices=['MAF', 'VCF', 'TCGI'], default='MAF')  # "Input format: MAF, VCF, or TCGI"
        optional_group.add_argument('--cohorts-file', type=str, help="Location of tar.gz container or directory for cohorts", default="cohorts.tar.gz", nargs='?')
        optional_group.add_argument('--cohort', "-c", type=str, help="Name of cohort with observed mutations", nargs='?', const="")

        advanced_group = parser.add_argument_group('Advanced arguments')
        advanced_group.add_argument('--profile', "-p", help="Override profile to calculate mutability, may also describe cohort size", type=str)
        advanced_group.add_argument('--nsamples', "-n", type=int, help="Override cohort size")
        advanced_group.add_argument('--threshold-driver', "-td", help="BScore threshold between Driver and Pontential Driver mutations", type=float, default=THRESHOLD_DRIVER)
        advanced_group.add_argument('--threshold-passenger', "-tp", help="BScore threshold between Pontential Driver and Passenger mutations", type=float, default=THRESHOLD_PASSENGER)

    def callback(self, args):
        """ Rank mutations """

        # no need to read mutations if there is a simple bail out
        if args.cohort is not None:
            if not Path(args.cohorts_file).is_file():
                logger.warning("Cohorts file missing. Download with \"mutagene fetch_cohorts\"")
                return
            if len(args.cohort) == 0:
                logger.warning("List of available cohorts:\n" + list_cohorts_in_tar(args.cohorts_file))
                return

        try:
            mutations, _, processing_stats = read_mutations(args.input_format, args.infile, args.genome, window_size=1)
        except Exception as e:
            e_message = getattr(e, 'message', repr(e))
            logger.warning(
                "Parsing {0} failed. "
                "Check that the input file is in {0} format "
                "or specify a different format using option -f \n"
                "{1}".format(args.input_format, e_message))

            if logger.root.level == logging.DEBUG:
                raise
            return

        msg = "Loaded {} mutations".format(processing_stats['loaded'])
        if processing_stats['skipped'] > 0:
            msg += " skipped {} mutations".format(processing_stats['skipped'])
        logger.info(msg)

        if not len(mutations) or not len(mutations[list(mutations.keys())[0]]):
            logger.warning('No mutations to rank. Check that the input file is in MAF format and correct genome assembly is chosen')
            return

        # read protein mutatations:
        args.infile.seek(0)
        protein_mutations, processing_stats = read_protein_mutations_MAF(args.infile, args.genome)
        msg = "Loaded {} protein mutations in {} samples".format(processing_stats['loaded'], processing_stats['nsamples'])
        if processing_stats['skipped'] > 0:
            msg += " skipped {} protein mutations".format(processing_stats['skipped'])
        logger.info(msg)

        cohort_size = len(mutations.keys())
        logger.info("Cohort size: {}".format(cohort_size))
        profile = get_pooled_multisample_mutational_profile(mutations, counts=True)
        cohort_aa_mutations = None

        # optionals:
        # overriding profile, cohort size and observed mutations based on precalc cohort
        if args.cohort is not None:
            profile, cohort_size, cohort_aa_mutations, cohort_na_mutations = read_cohort_mutations_from_tar(args.cohorts_file, args.cohort)
            logger.info("Precalculated cohort and profile loaded")
            logger.info("Cohort size: {}".format(cohort_size))

        # overriding profile:
        if args.profile:
            profile = read_profile_file(args.profile)
            if profile:
                logger.info('Profile overridden')
            else:
                logger.info('Profile could not be loaded')
                return
            cohort_size_new = read_cohort_size_from_profile_file(args.profile)
            if cohort_size_new:
                cohort_size = cohort_size_new
                logger.info('Cohort size loaded from profile N=' + str(cohort_size))

        # overriding n samples
        if args.nsamples:
            cohort_size = args.nsamples
            logger.info('Cohort size overridden N=' + str(cohort_size))

        logger.info("THRESHOLD_DRIVER: {}".format(args.threshold_driver))
        logger.info("THRESHOLD_PASSENGER: {}".format(args.threshold_passenger))

        # ranking:
        rank(protein_mutations, args.outfile, profile, cohort_aa_mutations, cohort_size, args.threshold_driver, args.threshold_passenger)
