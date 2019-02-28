# from collections import defaultdict
# import json
# import multiprocessing
# import requests

# import numpy as np
# import twobitreader as tbr

# from .dna import nucleotides, complementary_nucleotide
from .io import read_mutations, get_mutational_profile, write_profile_file


def calc_profile(infile, outfile, genome):
    mutations, processing_stats = read_mutations(infile, None, genome)
    profile = get_mutational_profile(mutations, counts=True)
    # print(profile)
    write_profile_file(outfile, profile)
    # print(profile)
