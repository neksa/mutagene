import re
import math
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import numpy as np
import pandas as pd

from tqdm import tqdm
# import functools
import pprint

from mutagene.dna import (
    nucleotides, complementary_nucleotide,
    bases_dict,
    # comp_dict,
    extended_nucleotides, complementary_extended_nucleotide)

from mutagene.io.motifs import get_known_motifs

import logging
logger = logging.getLogger(__name__)


def identify_motifs(samples_mutations, custom_motif=None, strand=None, threshold=None, dump_matches=None, stat_type=None):
    """
    :param samples_mutations: list of mutations from input file
    :param custom_motif: specified motif to search for
    :param strand: strand(s) to search on (T: transcribed, N: non-transcribed, A: any, or a combination theirof: 'TNA')
    :param dump_matches: pass through to process_mutations, stores all motif matches
    :param stat_type: pass through to process_mutations, choose statistical test
    :return: command-line output
    """
    motif_matches = []
    sig_motif_matches = []
    pvals = []

    if strand is None:
        strand = 'A'
    else:
        strand = set(strand)  # in case TNA codes repeat

    if threshold is None:
        threshold = 0.05

    if custom_motif:
        search_motifs = scanf_motif(custom_motif)
    else:
        motifs = get_known_motifs()
        search_motifs = motifs.copy()
    # search_motifs.extend(scanf_motif(custom_motif))

    _strand_map = {
        'T': 'transcribed',
        'N': 'non-transcribed',
        'A': 'any strand'
    }

    disable_progress_bar = logger.getEffectiveLevel() == logging.DEBUG

    for sample, mutations in tqdm(samples_mutations.items(), leave=False, disable=disable_progress_bar):
        if mutations is not None and len(mutations) > 0:
            first_mut_seq_with_coords = mutations[0][-1]
            window_size = (len(first_mut_seq_with_coords) - 1) // 2

            for m in tqdm(search_motifs, leave=False, disable=disable_progress_bar):
                for s in strand:
                    result, saved_data = process_mutations(
                        mutations,
                        m['motif'],
                        m['position'],
                        m['ref'],
                        m['alt'],
                        window_size,
                        s,
                        stat_type=stat_type)

                    if dump_matches:
                        for chrom, pos in saved_data['mutation_motif']:
                            dump_matches.write(
                                "chr{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    chrom, pos, int(pos) + 1, sample, m['logo'], _strand_map[s]))

                    debug_data = {
                        'sample': sample,
                        'motif': m['logo'],
                        'strand': s}
                    debug_data.update(result)
                    debug_string = pprint.pformat(debug_data, indent=4)
                    logger.debug(debug_string)

                    motif_matches.append({
                        'sample': sample,
                        'mutagen': m['name'],
                        'motif': m['logo'],
                        'strand': _strand_map[s],
                        'enrichment': result['enrichment'],
                        'mut_min': result['mutation_load'],
                        'mut_max': result['bases_mutated_in_motif'],
                        'odds_ratio': result['odds_ratio'],
                        'pvalue': result['pvalue']
                    })
                    pvals.append(result['pvalue'])

    qvalues = get_corrected_pvalues(pvals)
    for i, motif_dict in enumerate(motif_matches):
        motif_matches[i]['qvalue'] = qvalues[i]
        if motif_dict['mut_min'] == 0:
            continue
        if motif_dict['qvalue'] >= threshold:
            continue
        sig_motif_matches.append(motif_dict)

    return sig_motif_matches


def scanf_motif(custom_motif):
    """ recognize motif syntax like A[C>T]G and create a motif entry """
    m = re.search(
        r'([' + extended_nucleotides + ']*)\\[([' + nucleotides + '])>([' + extended_nucleotides + '])\\]([' + extended_nucleotides + ']*)',
        custom_motif.upper())
    if m:
        g = m.groups('')
        # print("GROUPS", m.group(1), m.group(2), m.group(3), m.group(4))
        entry = {}
        entry['logo'] = m.group(0)
        entry['motif'] = g[0] + g[1] + g[3]
        entry['position'] = len(g[0])
        entry['ref'] = g[1]
        entry['alt'] = g[2]
        if entry['ref'] == entry['alt']:
            return []
        entry['name'] = 'Custom motif'
        entry['references'] = ''
        return [entry, ]
    return []


def calculate_RR(ct):
    """
    Mutation is treatment
    No mutation is placebo

    :param ct: mutually exclusive counts of mutated matching motifs, matching mutations, matching motifs, and matching bases
    :return: enrichment or risk ratio
    """
    try:
        RR = ((ct.loc['mutation', 'motif'] / (ct.loc['mutation', 'motif'] + ct.loc['mutation', 'no motif'])) /
              (ct.loc['no mutation', 'motif'] / (ct.loc['no mutation', 'motif'] + ct.loc['no mutation', 'no motif'])))
    except ZeroDivisionError:
        RR = 0.0
    return RR


def calculate_RR_for_motif(ct):
    """
    Motif is treatment
    No motif is placebo

    :param ct: mutually exclusive counts of mutated matching motifs, matching mutations, matching motifs, and matching bases
    :return: enrichment or risk ratio
    """
    try:
        RR = ((ct.loc['mutation', 'motif'] / (ct.loc['mutation', 'motif'] + ct.loc['no mutation', 'motif'])) /
              (ct.loc['mutation', 'no motif'] / (ct.loc['mutation', 'no motif'] + ct.loc['no mutation', 'no motif'])))
    except ZeroDivisionError:
        RR = 0.0
    return RR


def calculate_OR(ct):
    """
    :param ct: mutually exclusive counts of mutated matching motifs, matching mutations, matching motifs, and matching bases
    :return: odds ratio
    """
    try:
        OR = (
            (ct.loc['mutation', 'motif'] / ct.loc['mutation', 'no motif']) /
            (ct.loc['no mutation', 'motif'] / ct.loc['no mutation', 'no motif']))
    except ZeroDivisionError:
        OR = 0.0
    return OR


def Haldane_correction(ct, pseudocount=0.5):
    """
    :param ct: mutually exclusive counts of mutated matching motifs, matching mutations, matching motifs, and matching bases
    :return: contigency table after Haldane correction is applied
    """
    """    apply Haldane correction (+ 0.5) if any of the values in the contingency table is zero """

    return ct + pseudocount if np.any(np.isclose(ct.to_numpy(), 0.0)) else ct


def calculate_mutation_load(N_mutations, enrichment):
    """
    Mutation load (minimum estimate) calculation following Gordenin et al protocol
    However, at this point motif matches are not filtered for p-value significance
    That's done in the end after multiple testing correction
    """
    mutation_load = 0.0
    if enrichment > 1.0:
        mutation_load = N_mutations * (enrichment - 1) / enrichment
    # elif p_value < p_value_threshold:   tests for enrichment depletion
    return mutation_load


def get_stats(ct, stat_type='fisher'):
    """
    Calculate Fisher and Chi2 test pvalues,
    :param ct: counts of mutated matching motifs, matching mutations, matching motifs, and matching bases
    :param stat_type: Type of pvalue (Fisher's ('fisher') or Chi-Square ('chi2'))
    :return: pvalue of the corresponding statistical test
    """
    p_val = 1.0

    if stat_type is None:
        stat_type = 'fisher'

    stat_type = stat_type.lower()

    acceptable_tests = ('fisher', 'chi2')
    if stat_type not in acceptable_tests:
        logger.warning('get_stats() can only calculate p-values for ' + str(acceptable_tests))

    if stat_type == 'fisher':
        try:
            p_val = stats.fisher_exact(ct, alternative="greater")[1]
            # if p_val > 0.05:
            #     p_val = stats.fisher_exact(ct, alternative="less")[1] #calculates if motif is underrepresented
        except ValueError:
            p_val = 1.0
    elif stat_type == 'chi2':
        try:
            p_val = stats.chi2_contingency(ct)[1]
        except ValueError:
            p_val = 1.0
    return p_val


def get_corrected_pvalues(p_values):
    qvalues = []
    if len(p_values):
        qvalues = multipletests(pvals=p_values, method='fdr_bh')[1]
    return qvalues


# @functools.lru_cache(maxsize=None)
def get_rev_comp_seq(sequence):
    """
    :param sequence: forward DNA sequence
    :return: reverse complimentary DNA sequence
    """
    # rev_comp_seq = "".join([complementary_nucleotide[i] for i in reversed(sequence)])
    cn = complementary_nucleotide
    return [(i[0], i[1], cn[i[2]], '-') for i in reversed(sequence)]


def mutated_base(mutation, ref, alt):
    """
    :param mutation: [(record.CHROM, record.POS, record.REF, record.ALT)]
    :param ref: list the nucleotide base pre-mutation
    :param alt: list the nucleotide base post-mutation
    :return: True if mutation matches the specified ref and alt
    """
    # makes sure single base substitution
    _, _, mut_ref, mut_alt = mutation
    if mut_alt and mut_ref and len(mut_ref) == 1 and len(mut_alt) == 1 and mut_ref != mut_alt:
        # mutation matches the substitution
        if mutation[2] in bases_dict[ref] and mutation[3] in bases_dict[alt]:
            return True


def find_matching_motifs(seq, motif, motif_position):
    """
    :param seq: DNA sequence
    :param motif: specified motif
    :param motif_position: position of mutated base in motif, 0-base numbering
    :return: generator of matching positions

    TODO: SLOW algorithm O(n * m). Need to create a suffix tree with regexp
    """
    # print("Looking for motif {} in {}, {}".format(motif, sequence, len(sequence) - len(motif)))
    for i in range(len(seq) - len(motif) + 1):
        # s = seq[i: i + len(motif)]
        # print(s)
        for j, c in enumerate(motif):
            if seq[i + j][2] not in bases_dict[c]:
                break
        else:
            yield seq[i + motif_position]


def find_matching_bases(seq, ref, motif, motif_position):
    """
    :param seq:
    :param ref:
    :param motif:
    :param motif_position:
    :return: bases that match mutations
    """
    for i in range(motif_position, len(seq) - (len(motif) - motif_position) + 1):
        # range excludes border of sequence that may be motifs that don't fit window size
        if seq[i][2] in bases_dict[ref]:
            yield seq[i]


def make_contingency_table(
        array=None,
        motif_mutation=None,
        no_motif_mutation=None,
        motif_no_mutation=None,
        no_motif_no_mutation=None):
    """ Make a 2x2 contingency table out of a numpy array or four integers"""

    if array is not None:
        assert isinstance(array, np.ndarray)
        assert array.shape == (2, 2)
    else:
        array = np.array([
            [motif_mutation, no_motif_mutation],
            [motif_no_mutation, no_motif_no_mutation]
        ])
    contingency_table = pd.DataFrame(array)
    contingency_table.columns = ["motif", "no motif"]
    contingency_table.index = ["mutation", "no mutation"]
    return contingency_table


def process_mutations(mutations, motif, motif_position, ref, alt, range_size, strand, stat_type=None):
    """
    :param mutations: mutations to be analyzed
    :param motif: specified motif to search for
    :param motif_position: location of mutation in motif, 0-base numbering from left of motif
    :param ref: base pre-mutation
    :param alt: base post-mutation
    :param range_size: how far in the motif to search for
    :param strand: strand motif should be searched on
    :param stat_type: type of pvalue: Fisher's (default) or Chi-Square
    :param dump_matches: an optional file handle to save all mutations matching the motif regardless of their significance
    :return: (results summary disctionary, data_dump with stored_data or None if dump_matches is None)
    """
    assert range_size >= 0
    assert len(ref) == 1
    assert len(alt) == 1
    assert 0 <= motif_position < len(motif)
    assert len(set(strand) - set("ATN")) == 0, "[process_mutations] only A, T, N allowed in strand parameter"

    matching_bases = set()
    matching_motifs = set()
    matching_mutated_motifs = set()
    matching_mutated_bases = set()

    # extra loop for sample in sample list
    for chrom, pos, transcript_strand, x, y, seq in mutations:
        # extract the longest sequence we would ever need (motif + range_size); range size = # bases outside mutation
        mutation = chrom, pos, x, y
        rev_seq = get_rev_comp_seq(seq)

        # assuming that all mutations are reported in '+' reference strand
        if strand == 'A' or (strand == 'T' and transcript_strand == '+') or (strand == 'N' and transcript_strand == '-'):
            # not mutated:
            for ref_match in find_matching_bases(seq, ref, motif, motif_position):
                matching_bases.add(ref_match[0:2])
            for motif_match in find_matching_motifs(seq, motif, motif_position):
                matching_motifs.add(motif_match[0:2])

            # mutated:
            if mutated_base(mutation, ref, alt):
                # m = (mutation[0], mutation[1], mutation[2], "+")
                matching_mutated_bases.add(mutation[0:2])

                context_of_mutation = seq[range_size - motif_position: range_size - motif_position + len(motif)]
                for motif_match in find_matching_motifs(context_of_mutation, motif, motif_position):
                    matching_mutated_motifs.add(motif_match[0:2])

        if strand == 'A' or (strand == 'T' and transcript_strand == '-') or (strand == 'N' and transcript_strand == '+'):
            # rev compl: not mutated:
            for ref_match in find_matching_bases(rev_seq, ref, motif, motif_position):
                matching_bases.add(ref_match[0:2])

            for motif_match in find_matching_motifs(rev_seq, motif, motif_position):
                matching_motifs.add(motif_match[0:2])

            # rev compl: mutated:
            if mutated_base(mutation, complementary_extended_nucleotide[ref], complementary_extended_nucleotide[alt]):
                # m = (mutation[0], mutation[1], mutation[2], "-")
                matching_mutated_bases.add(mutation[0:2])

                # rev comp:
                context_of_mutation = rev_seq[range_size - motif_position: range_size - motif_position + len(motif)]
                for motif_match in find_matching_motifs(context_of_mutation, motif, motif_position):
                    matching_mutated_motifs.add(motif_match[0:2])

    motif_mutation_count = len(matching_mutated_motifs)  # bases mutated in motif
    stat_mutation_count = len(matching_mutated_bases - matching_mutated_motifs)  # bases mutated not in motif
    stat_motif_count = len(matching_motifs - matching_mutated_motifs)  # bases not mutated in motif
    stat_ref_count = len(matching_bases - (matching_motifs | matching_mutated_bases))  # bases not mutated not in motif

    # number of A[T>G]T occurrences               motif_mutation_count
    # / number of [T>G] occurrences               stat_mutation_count + motif_mutation_count
    # ----------
    # number of ATT occurrences in DNA context    stat_motif_count
    # / number of T occurrences in DNA context    stat_ref_count + stat_motif_count

    contingency_table = make_contingency_table(
        motif_mutation=motif_mutation_count,
        no_motif_mutation=stat_mutation_count,
        motif_no_mutation=stat_motif_count,
        no_motif_no_mutation=stat_ref_count)

    # data={
    # "'{}>{}' mutation".format(ref, alt): [stat_mutation_count, motif_mutation_count],
    # "no '{}>{}' mutation".format(ref, alt): [stat_ref_count, stat_motif_count]},
    # index=("no '{}' motif".format(motif), "'{}' motif".format(motif)))
    logger.debug("\n" + contingency_table.to_string() + "\n")

    logger.debug("({} / ({} + {})  ) /  ({} / ({} + {}))".format(
        contingency_table.loc['mutation', 'motif'],
        contingency_table.loc['mutation', 'motif'],
        contingency_table.loc['mutation', 'no motif'],
        contingency_table.loc['no mutation', 'motif'],
        contingency_table.loc['no mutation', 'motif'],
        contingency_table.loc['no mutation', 'no motif']))

    contingency_table = Haldane_correction(contingency_table)

    enrichment = risk_ratio = calculate_RR(contingency_table)  # enrichment = risk ratio
    odds_ratio = calculate_OR(contingency_table)

    p_val = get_stats(contingency_table, stat_type)

    mut_load = calculate_mutation_load(motif_mutation_count, enrichment)

    result = {
        'enrichment': enrichment,  # AKA risk ratio
        'odds_ratio': odds_ratio,
        'mutation_load': math.ceil(mut_load),
        'pvalue': p_val,
        'bases_mutated_in_motif': motif_mutation_count,
        'bases_mutated_not_in_motif': stat_mutation_count,
        'bases_not_mutated_in_motif': stat_motif_count,
        'bases_not_mutated_not_in_motif': stat_ref_count,
        'total_mutations': len(mutations)
    }

    saved_matches = {
        'mutation_motif': matching_mutated_motifs
    }
    return result, saved_matches
