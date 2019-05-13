#!/usr/bin/env python
# -*- coding: utf-8 -*- 

# from collections import defaultdict
# from pprint import pprint
# from statsmodels.stats.multitest import multipletests

import scipy.stats as stats
import math
import numpy as np
# from itertools import cycle

from mutations_io import *
from .annotated_motifs import annotated_motifs


def identify_motifs(mutations, custom_motifs=None):
    motif_matches = []

    selected_motifs = custom_motifs if custom_motifs else motifs

    if mutations is not None and len(mutations) > 0:
        first_mut_seq_with_coords = mutations[0][-1]
        window_size = (len(first_mut_seq_with_coords) - 1) // 2

        for m in selected_motifs:
            print("IDENTIFYING MOTIF: ", m['name'])
            result = get_enrichment(mutations, m['motif'], m['position'], m['ref'], m['alt'], window_size)
            print(result)
            print()
            # enrichment, math.ceil(mut_load), p_val[1], p_val[2], motif_mutation_count, mutation_count, motif_count, ref_count
            enrichment, mut_burden, pvalue_fisher, pvalue_chi, motif_mut, mut, motif_count, ref_count = result

            motif_matches.append({
                'name': m['name'],
                'motif': m['logo'],
                'enrichment': enrichment,
                'pvalue': pvalue_fisher,
                'mutations': mut_burden,
                })
    return motif_matches


def get_stats(motif_mutation_count, mutation_count, motif_count, ref_count):
    p_val = stats.fisher_exact([[mutation_count, motif_mutation_count], [ref_count, motif_count]], alternative= "less")
    p_value = p_val[1]
    chi_array = np.array([[mutation_count , motif_mutation_count], [ref_count, motif_count]])
    if motif_mutation_count > 0 and mutation_count > 0 and motif_count > 0 and ref_count > 0:
        chi = stats.chi2_contingency(chi_array)[1]
    else:
        chi = 1.00

    # if p_value <= 0.05:
    #  qvalues = multipletests(pvals=p_value, method='fdr_bh')
    #  print(qvalues)
    #  if qvalues[3] <= 0.05:
    #     print("significant")
    # print("odds_ratio: ", "p-value")
    return p_val[0], p_value, chi


def get_rev_comp_seq(sequence):
    # rev_comp_seq = "".join([complementary_nucleotide[i] for i in reversed(sequence)])
    rev_comp_seq = [(i[0], i[1], complementary_nucleotide[i[2]], "-") for i in reversed(sequence)]
    return rev_comp_seq


def mutated_base(mutation, ref, alt):
    """
    :param mutation: [(record.CHROM, record.POS, record.REF, record.ALT)]
    :param ref: the nucleotide base pre-mutation
    :param alt: the nucleotide base post-mutation
    :return: how many mutations match specified ref and alt
    """
    assert ref != alt, "mutation should have different ref and alt nucleotides"

    # makes sure single base substitution
    if mutation[3] and mutation[2] and len(mutation[2]) == 1 \
            and len(mutation[3]) == 1 and len(mutation[3]) == 1 \
            and len(mutation[2]) == 1 and mutation[2] != mutation[3]:
        # mutation matches the substitution
        if mutation[2] in ref and mutation[3] in alt:
            return True


def find_matching_motifs(seq, motif, motif_position):
    # print("Looking for motif {} in {}, {}".format(motif, sequence, len(sequence) - len(motif)))
    for i in range(len(seq) - len(motif) + 1):
        s = seq[i: i + len(motif)]
        for j, char in enumerate(motif):
            if s[j][2] not in bases_dict[char]:
                break
        else:
            yield seq[i + motif_position]


def find_matching_bases(seq, ref, motif, motif_position):
    # print("Looking for motif {} in {}, {}".format(motif, sequence, len(sequence) - len(motif)))
    for i in range(len(seq)):
        s = seq[i]
        if s[2] in ref:
            yield seq[i]


def get_enrichment(mutations, motif, motif_position, ref, alt, range_size):
    assert range_size >= 0
    assert range_size == 50, "code will break if you use different range size"
    assert len(ref) == 1
    assert len(alt) == 1
    assert 0 <= motif_position < len(motif)

    matching_bases = set()
    matching_motifs = set()
    matching_mutated_motifs = set()
    matching_mutated_bases = set()
    matching_mutated_refs = set()

    # extra loop for sample in sample list
    for chrom, pos, x, y, seq in mutations:
        # extract the longest sequence we would ever need (motif + range_size)
        mutation = chrom, pos, x, y
        rev_seq = get_rev_comp_seq(seq)

        rev_motif_pos = len(motif) - motif_position - 1
        ref_seq = seq[motif_position:len(seq) - (len(motif) - motif_position) + 1]

        if mutation[2] == bases_dict[ref] or mutation[2] == comp_dict[ref]:  # checks to see how many times ref base is mutated
            matching_mutated_refs.add(mutation[0:2])

        # not mutated:
        for ref_match in find_matching_bases(ref_seq, ref, motif, motif_position):
            matching_bases.add(ref_match[0:2])

        for motif_match in find_matching_motifs(seq, motif, motif_position):
            matching_motifs.add(motif_match[0:2])

        # rev compl: not mutated:
        ref_rev_seq = rev_seq[rev_motif_pos: len(seq) - (len(motif) - rev_motif_pos) + 1]

        for ref_match in find_matching_bases(ref_rev_seq, ref, motif, len(motif) - motif_position - 1):
            matching_bases.add(ref_match[0:2])

        for motif_match in find_matching_motifs(rev_seq, motif, len(motif) - motif_position - 1):
            matching_motifs.add(motif_match[0:2])

        # mutated:
        if mutated_base(mutation, bases_dict[ref], bases_dict[alt]):
            # m = (mutation[0], mutation[1], mutation[2], "+")
            matching_mutated_bases.add(mutation[0:2])

            context_of_mutation = seq[range_size - motif_position: range_size - motif_position + len(motif)]
            for motif_match in find_matching_motifs(context_of_mutation, motif, motif_position):
                matching_mutated_motifs.add(motif_match[0:2])

        if mutated_base(mutation, comp_dict[ref], comp_dict[alt]):
            # m = (mutation[0], mutation[1], mutation[2], "-")
            matching_mutated_bases.add(mutation[0:2])

            # rev comp:
            context_of_mutation = rev_seq[range_size - motif_position: range_size - motif_position + len(motif)]
            for motif_match in find_matching_motifs(context_of_mutation, motif, len(motif) - motif_position - 1):
                matching_mutated_motifs.add(motif_match[0:2])

    motif_mutation_count = len(matching_mutated_motifs)  # mutated
    mutation_count = len(matching_mutated_bases - matching_mutated_motifs)  # mutated
    ref_count = len(matching_bases - matching_motifs - matching_mutated_refs)  # not mutated
    motif_count = len(matching_motifs - matching_mutated_refs)  # not mutated

    stat_motif_count = len(matching_motifs - matching_mutated_motifs - matching_mutated_refs)
    stat_ref_count = len(matching_bases - matching_motifs - matching_mutated_bases - matching_mutated_refs)

    try:
        enrichment = (float(motif_mutation_count) / float(mutation_count)) / (float(motif_count) / float(ref_count))
    except ZeroDivisionError:
        enrichment = 0.0
    p_val = get_stats(motif_mutation_count, mutation_count, stat_motif_count, stat_ref_count)
    chi = p_val[2]

    if enrichment > 1 and p_val[1] < 0.05 and chi < 0.05:
        mut_load = (motif_mutation_count * (enrichment - 1)) / enrichment
    else:
        mut_load = 0.0
    if chi == 1.00:
        chi = "chi cannot compute due to zero count"
    return enrichment, math.ceil(mut_load), p_val[1], chi, motif_mutation_count, mutation_count, motif_count, ref_count


def scan_window_range(filename, motif, motif_position, ref, alt, assembly=None):
    mut_burden_count = []
    range_sizes = []
    with open("APOBEC_window.csv", "w") as csvfile: #create csv file
        filewriter = csv.writer(csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow([motif])
        filewriter.writerow(["range_size", "enrichment", "motif_mut_count", "mut_count", "motif_count", "ref_count", "mutation load", "fisher p-value", "chi p-value"])
        for val in range(0, 200 + 1, 5):
            range_sizes.append(val * 2 + len(motif))
            analysis = get_enrichment(filename, motif, motif_position, ref, alt, val, assembly=assembly)
            result = analysis[0]
            mut_burden = analysis[1]
            mut_burden_count.append(mut_burden)
            pvalue = analysis[2]
            chi_pvalue = analysis[3]
            motif_mut = analysis[4]
            mut = analysis[5]
            motif_count = analysis[6]
            ref_count = analysis[7]
            filewriter.writerow([val, result, motif_mut, mut, motif_count, ref_count, mut_burden, pvalue, chi_pvalue])
        plt.plot(range_sizes, mut_burden_count, color='b')
        plt.xlabel("Window Size")
        plt.ylabel("Mutational Burden")
        plt.title("APOBEC Mutational Burden Depends Upon Window Size")
        plt.show()
        return list(zip(range_sizes, mut_burden_count))


if __name__ == '__main__':
    pass
    # pprint(get_values("tcga_A0C8.maf", "TCW", 1, "C", "G", assembly=37))
    # mutations = read_MAF(filename)
    # pprint(get_enrichment(mutations, "TW", 0, "T", "G", 29, assembly=37))
    # pprint(get_enrichment("tcga_A0C8.maf", "NTW", 1, "T", "G", 29, assembly=37))

    # pprint(get_enrichment("tcga_A0C8.maf", "NTW", 1, "T", "G", 29, assembly=37))
    # pprint(get_enrichment("tcga_A0C8.maf", "TCW", 1, "T", "G", 29, assembly=37))
    # pprint(get_enrichment("tcga_A0C8.maf", "AGA", 1, "G", "A", 29, assembly=37))
    # pprint(get_enrichment("tcga_A0C8.maf", "TGC", 1, "G", "T", 29, assembly=37))
