# from collections import defaultdict
# from mutagene.dna import complementary_nucleotide

import logging
from collections import Counter
from functools import reduce
from operator import add

import numpy as np
from numpy.random import multinomial

from mutagene.dna import complementary_nucleotide
from mutagene.io.context_stats import (
    merge_context_stats,
    new_context_stats,
    report_assembly_mismatch,
)
from mutagene.io.context_window import read_mutations
from mutagene.io.matrix import channel_labels
from mutagene.io.mutations_profile import read_auto_profile
from mutagene.io.profile import get_profile_attributes_dict, write_profile_file

logger = logging.getLogger(__name__)


def calc_profile(infile, outfile, genome, fmt="auto", keep_filtered=False, filter_column="FILTER"):
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
            f, fmt=fmt, asm=genome, keep_filtered=keep_filtered, filter_column=filter_column
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
    if not counts and total_mut_number == 0:
        # Nothing to take a share of. Dividing produced ZeroDivisionError from
        # inside the loop, several frames from the empty input that caused it.
        raise ValueError("Cannot compute a frequency profile from zero mutations")
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


def read_sample_profiles(fmt, stream, genome, context, *, keep_filtered, filter_column):
    """Count substitutions separately for each sample in the requested sequence context."""
    if context not in {"SBS96", "SBS1536"}:
        raise ValueError("Substitution calls require an SBS context")
    mutations, sites, _ = read_mutations(
        fmt,
        stream,
        genome,
        window_size=2 if context == "SBS1536" else 1,
        keep_filtered=keep_filtered,
        filter_column=filter_column,
    )
    if context == "SBS96":
        return get_multisample_mutational_profile(mutations, counts=True)

    labels = channel_labels(context)
    profiles = {}
    for sample, records in sites.items():
        counts = Counter()
        for _, _, _, ref, alt, window in records:
            sequence = "".join(base for _, _, base, _ in window)
            if len(sequence) != 5 or any(base not in "ACGT" for base in sequence):
                continue
            if sequence[2] != ref:
                ref, alt = complementary_nucleotide[ref], complementary_nucleotide[alt]
            if sequence[2] != ref or ref == alt:
                continue
            if ref in "AG":
                sequence = "".join(complementary_nucleotide[base] for base in reversed(sequence))
                alt = complementary_nucleotide[alt]
            counts[f"{sequence[:2]}[{sequence[2]}>{alt}]{sequence[3:]}"] += 1
        if not counts:
            raise ValueError(f"No substitutions have complete pentanucleotide contexts in {sample}")
        skipped = len(records) - sum(counts.values())
        if skipped:
            logger.warning(
                "Skip %d substitutions without a usable pentanucleotide in %s", skipped, sample
            )
        profiles[sample] = np.array([counts[label] for label in labels])
    return profiles


def get_pooled_multisample_mutational_profile(samples_mutations, counts=False):
    mutational_profile_dict = reduce(
        add, (Counter(dict(mutations)) for mutations in samples_mutations.values())
    )
    return get_mutational_profile(mutational_profile_dict, counts)


def generate_resampled_profiles(profile, k, rng=None):
    profile = np.array(profile)
    N = np.sum(profile)
    sampler = multinomial if rng is None else rng.multinomial
    new_profiles = sampler(int(N), profile / N, size=k)
    for i in range(k):
        yield new_profiles[i]
