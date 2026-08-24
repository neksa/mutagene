# from collections import defaultdict
# from mutagene.dna import complementary_nucleotide

import logging
from collections import Counter
from functools import reduce
from operator import add

import numpy as np
from numpy.random import multinomial

from mutagene.io.context_stats import (
    merge_context_stats,
    new_context_stats,
    report_assembly_mismatch,
)
from mutagene.io.mutations_profile import read_auto_profile
from mutagene.io.profile import get_profile_attributes_dict, write_profile_file

logger = logging.getLogger(__name__)


def calc_profile(infile, outfile, genome, fmt="auto", keep_filtered=False):
    """Write a 96-channel profile for the given inputs.

    Returns the combined processing stats, including the ``ref_mismatches``
    count that identifies a wrong genome assembly. Callers that only want the
    file may ignore it.
    """
    all_mutations = {}
    combined_stats = new_context_stats()
    combined_stats.update({"loaded": 0, "skipped": 0})
    for f in infile:
        mutations, processing_stats = read_auto_profile(
            f, fmt=fmt, asm=genome, keep_filtered=keep_filtered
        )
        msg = "Loaded {} mutations".format(processing_stats["loaded"])
        if processing_stats["skipped"] > 0:
            msg += " skipped {} mutations due to mismatches with the reference genome".format(
                processing_stats["skipped"]
            )
        logger.info(msg)
        merge_context_stats(combined_stats, processing_stats)
        combined_stats["loaded"] += processing_stats.get("loaded", 0)
        combined_stats["skipped"] += processing_stats.get("skipped", 0)
        all_mutations = {
            k: all_mutations.get(k, 0) + mutations.get(k, 0)
            for k in set(all_mutations) | set(mutations)
        }

    report_assembly_mismatch(combined_stats, combined_stats["loaded"], genome)

    if sum(all_mutations.values()) == 0:
        logger.warning("Can not create profile")
        return combined_stats
    profile = get_mutational_profile(all_mutations, counts=True)
    write_profile_file(outfile, profile)
    return combined_stats


def get_mutational_profile(mutational_profile_dict, counts=False):
    attrib = get_profile_attributes_dict()
    values = []
    total_mut_number = sum(mutational_profile_dict.values())
    for i, attr in enumerate(attrib):
        number = mutational_profile_dict.get(attr["context"] + attr["mutation"], 0)
        # freq = 0.000001 * number / total_mut_number
        if counts:
            freq = number
        else:
            freq = number / float(total_mut_number)
        # trinucleotide = attr['context'][0] + attr['mutation'][0] + attr['context'][1]
        # trinucleotide_freq = exome_trinucleotide_freq[trinucleotide]
        # values.append(3.0 * freq / trinucleotide_freq)
        values.append(freq)
    return values


def get_multisample_mutational_profile(samples_mutations, counts=False):
    samples_profiles = {}

    for sample, mutations in samples_mutations.items():
        samples_profiles[sample] = get_mutational_profile(mutations, counts)

    return samples_profiles


def get_pooled_multisample_mutational_profile(samples_mutations, counts=False):
    mutational_profile_dict = reduce(
        add, (Counter(dict(mutations)) for mutations in samples_mutations.values())
    )
    return get_mutational_profile(mutational_profile_dict, counts)


def generate_resampled_profiles(profile, k):
    profile = np.array(profile)
    N = np.sum(profile)
    new_profiles = multinomial(N, profile / N, size=k)
    for i in range(k):
        yield new_profiles[i]
