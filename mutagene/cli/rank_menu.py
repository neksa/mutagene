import argparse
import sys
import logging
from pathlib import Path

from mutagene.mutability.mutability import rank, THRESHOLD_DRIVER, THRESHOLD_PASSENGER
from mutagene.io.cohorts import read_cohort_mutations_from_tar
from mutagene.io.cohorts import read_cohort_size_from_profile_file, list_cohorts_in_tar
from mutagene.io.profile import read_profile_file
from mutagene.io.protein_mutations_MAF import read_MAF_with_genomic_context


logger = logging.getLogger(__name__)
genome_error_message = """requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more
                          Use mutagene fetch to download genome assemblies"""


class RankMenu(object):
    def __init__(self, parser):
        parser.add_argument("--infile", "-i", nargs='*', help="Input file in MAF format", type=argparse.FileType('r'))
        parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
        parser.add_argument('--genome', "-g", help="Location of genome assembly file", type=str)

        # parser.add_argument('--mode', "-m", help="n (nucleotide) or aa (amino acid) mode", type=str, default="aa")

        parser.add_argument('--cohorts-file', type=str, help="Location of tar.gz container or directory for cohorts", default="cohorts.tar.gz")
        parser.add_argument('--cohort', "-c", type=str, help="Name of cohort with observed mutations")

        parser.add_argument('--profile', "-p", help="Override profile to calculate mutability, may also describe cohort size", type=str)
        parser.add_argument('--nsamples', "-n", type=int, help="Override cohort size")
        parser.add_argument('--threshold-driver', "-td", help="BScore threshold between Driver and Pontential Driver mutations", type=float, default=THRESHOLD_DRIVER)
        parser.add_argument('--threshold-passenger', "-tp", help="BScore threshold between Pontential Driver and Passenger mutations", type=float, default=THRESHOLD_PASSENGER)

    @classmethod
    def callback(cls, args):
        # print('RankMenu', args.infile)
        # if not args.resource:
        #     logger.warning('No resource specified')
        #     sys.exit(1)

        # getattr(cls, args.resource)(args)

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
            if len(args.infile) > 1:
                logger.info('Multiple input files provided')
            mutations_to_rank = {}
            processing_stats = {'loaded': 0, 'skipped': 0}
            for infile in args.infile:
                mut, stats = read_MAF_with_genomic_context(infile, args.genome)
                mutations_to_rank.update(mut)
                processing_stats['loaded'] += stats['loaded']
                processing_stats['skipped'] += stats['skipped']
        else:
            mutations_to_rank, processing_stats = read_MAF_with_genomic_context(args.infile, args.genome)

        if not len(mutations_to_rank):
            logger.warning('MutaGene rank failed: No mutations to rank. Check that the infile is in MAF format')
            return

        msg = "Loaded {} mutations".format(processing_stats['loaded'])
        if processing_stats['skipped'] > 0:
            msg += " skipped {} mutations".format(processing_stats['skipped'])
        logger.info(msg)

        # mutations_to_rank = list(mutations_to_rank)

        # not_unique = len(mutations_to_rank)
        # mutations_to_rank = list(set(mutations_to_rank))
        # unique = len(mutations_to_rank)
        # if not_unique != unique:
        #     logger.info("Removed {} duplicate mutations".format(not_unique - unique))

        logger.info("THRESHOLD_DRIVER: {}".format(args.threshold_driver))
        logger.info("THRESHOLD_PASSENGER: {}".format(args.threshold_passenger))

        rank(mutations_to_rank, args.outfile, profile, cohort_aa_mutations, cohort_size, args.threshold_driver, args.threshold_passenger)


"""

# parser.add_argument('integers', metavar='N', type=int, nargs='+',
#                     help='an integer for the accumulator')
# parser.add_argument('--sum', dest='accumulate', action='store_const',
#                     const=sum, default=max,
#                     help='sum the integers (default: find the max)')
parser.add_argument("--infile", "-i", nargs='*', help="Input file format", type=argparse.FileType('r'))
parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
# parser.add_argument('--signatures', type=str)
# parser.add_argument('--motifs', type=str)

# parser.add_argument('--threshold', "-t", help="B-score thresholds for drivers", type=float)
# parser.add_argument('--mut-rate', "-r", help="Mutation rate overrides the rate inferred from profile", type=float)

"""
