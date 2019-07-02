
from mutagene.io.profile import write_profile_file, get_profile_attributes_dict
from mutagene.io.mutations_profile import read_auto_profile

import logging
logger = logging.getLogger(__name__)


def calc_profile(infile, outfile, genome):
    all_mutations = {}
    for f in infile:
        mutations, processing_stats = read_auto_profile(f, fmt='auto', asm=genome)
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


def get_mutational_profile(mutations, counts=False):
    attrib = get_profile_attributes_dict()
    values = []
    total_mut_number = sum(mutations.values())
    for i, attr in enumerate(attrib):
        number = mutations.get(attr['context'] + attr['mutation'], 0)
        # freq = 0.000001 * number / total_mut_number
        if counts:
            freq = number
        else:
            freq = number / float(total_mut_number)
        # trinucleotide = attr['context'][0] + attr['mutation'][0] + attr['context'][1]
        # trinucleotide_freq = exome_trinucleotide_freq[trinucleotide]
        # values.append(3.0 * freq / trinucleotide_freq)
        values.append(freq)
    # print(values)
    return values
