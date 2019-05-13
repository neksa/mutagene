# from collections import defaultdict
# import json
# import multiprocessing
# import requests

# import numpy as np
# import twobitreader as tbr

# from .dna import nucleotides, complementary_nucleotide
from mutagene.io.io import read_mutations, get_mutational_profile, write_profile_file

import logging
logger = logging.getLogger(__name__)


def calc_profile(infile, outfile, genome):
    all_mutations = {}
    for f in infile:
        mutations, processing_stats = read_mutations(f, None, genome)
        msg = "Loaded {} mutations".format(processing_stats['loaded'])
        if processing_stats['skipped'] > 0:
            msg += " skipped {} mutations due to mismatches with the reference genome".format(processing_stats['skipped'])
        logger.info(msg)
        all_mutations = {k: all_mutations.get(k, 0) + mutations.get(k, 0) for k in set(all_mutations) | set(mutations)}
    if sum(all_mutations.values()) == 0:
        logger.warn('Can not create profile')
        return
    profile = get_mutational_profile(all_mutations, counts=True)
    # print(profile)
    write_profile_file(outfile, profile)
    # print(profile)
