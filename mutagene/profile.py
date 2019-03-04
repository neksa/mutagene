# from collections import defaultdict
# import json
# import multiprocessing
# import requests

# import numpy as np
# import twobitreader as tbr

# from .dna import nucleotides, complementary_nucleotide
from .io import read_mutations, get_mutational_profile, write_profile_file

import logging
logger = logging.getLogger(__name__)


def calc_profile(infile, outfile, genome):
    mutations, processing_stats = read_mutations(infile, None, genome)
    msg = "Loaded {} mutations".format(processing_stats['loaded'])
    if processing_stats['skipped'] > 0:
        msg += " skipped {} mutations due to mismatches with the reference genome".format(processing_stats['skipped'])
    logger.info(msg)
    profile = get_mutational_profile(mutations, counts=True)
    # print(profile)
    write_profile_file(outfile, profile)
    # print(profile)
