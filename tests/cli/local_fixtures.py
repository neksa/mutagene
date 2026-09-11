"""Fixture paths that need no network.

Most CLI tests only need "a genome and a MAF that matches it". They used a
genome downloaded from UCSC, so a runner that could not reach it ran fewer
tests and reported lower coverage for the same commit, with nothing to say so.

sample_small.maf was generated from LOCAL_GENOME, so every reference allele
matches and no mutation is discarded.
"""

import os

LOCAL_GENOME = os.path.join(os.path.dirname(__file__), "..", "motifs", "data", "test_genome.2bit")
LOCAL_MAF = os.path.join(os.path.dirname(__file__), "data", "sample_small.maf")

# Mutations in sample_small.maf, all of which load against LOCAL_GENOME.
LOCAL_MAF_MUTATIONS = 24
LOCAL_MAF_SAMPLES = 3
