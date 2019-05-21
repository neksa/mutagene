import pandas as pd

from collections import defaultdict, OrderedDict
from functools import lru_cache, reduce

from scipy.stats import binom_test
from statsmodels.stats.multitest import multipletests

from mutagene.dna import *
from mutagene.io.profile import get_profile_attributes_dict

import logging
logger = logging.getLogger(__name__)


THRESHOLD_DRIVER = 8.030647e-05  # Max MCC
THRESHOLD_PASSENGER = 0.003440945  # 10% FPR


def get_all_codon_substitutions(P, seq5, Q):
    """
        Find all possible single nucleotide substitutions that mutate seq5 codon + context to amino acid Q
    """
    seq5_codon = seq5[1:-1]
    mutated_seq5s = []
    if seq5_codon not in codon_table:
        logger.warning("Invalid codon " + seq5_codon)
    if codon_table[seq5_codon] != P:
        logger.warning("Codon P does not match the sequence " + seq5_codon + " " + P)
    for codon, A in codon_table.items():
        if A != Q:
            continue
        n_changes = 0
        for i in range(3):
            if codon[i] != seq5_codon[i]:
                n_changes += 1
        if n_changes == 1:
            mutated_seq5s.append(seq5[0] + codon + seq5[-1])
    return mutated_seq5s


@lru_cache(maxsize=1000)
def predict_driver(observed, N, p, threshold_driver, threshold_passenger):
    """
    Binomial P-value based on the number of observed mutations, number of samples and expected mutability
    """
    # TODO: thresholds should be defined in the config
    # need more tumor samples than mutations for sure

    if N > 0 and observed > 0 and observed <= N:
        pvalue = binom_test(observed, N, p, alternative='greater')
    else:
        pvalue = 1.0

    if pvalue <= threshold_driver:
        prediction = "Driver"
        prediction_value = 1
    elif pvalue <= threshold_passenger:
        prediction = "Potential driver"
        prediction_value = 0.5
    elif pvalue < 1.0:
        prediction = "Passenger"
        prediction_value = 0
    else:
        prediction = "Undefined"
        prediction_value = -1

    return pvalue, prediction, prediction_value


def calculate_codon_mutability(mutation_model, seq5, mutated_seq5s):
    positional_mutability = defaultdict(float)
    for mutated_seq5 in mutated_seq5s:
        j = 0
        y = None
        for i in range(3):
            j = i + 1
            y = mutated_seq5[j]
            if seq5[j] != y:
                break
        if j > 0 and y is not None:
            # should always be > 0
            tri = seq5[j - 1: j + 2]
            positional_mutability[i] += mutation_model[tri][y]
    return 1.0 - reduce(lambda x, y: x * (1.0 - y), positional_mutability.values(), 1.0)


def calculate_base_substitution_mutability(counts_profile, cohort_size):
    assert cohort_size > 0
    assert len(counts_profile) == 96

    mutability = defaultdict(dict)
    counts_profile_dict = defaultdict(dict)

    for i, attrib in enumerate(get_profile_attributes_dict()):
        x, y = attrib['mutation']
        p5, p3 = attrib['context']
        tri = p5 + x + p3
        counts_profile_dict[tri][y] = counts_profile[i]

    for tri in counts_profile_dict.keys():
        tri_sum = float(sum(counts_profile_dict[tri].values()))
        mutability_tri = tri_sum / (float(cohort_size) * float(exome_trinucleotides[tri]))
        for y, y_counts in counts_profile_dict[tri].items():
            f = 0.0
            if tri_sum > 0:
                f = float(y_counts) / tri_sum
            mutability[tri][y] = f * mutability_tri
            mutability[complementary_trinucleotide[tri]][complementary_nucleotide[y]] = mutability[tri][y]

    return mutability


# import pprint
def rank(mutations_to_rank, outfile, profile, cohort_aa_mutations, cohort_size, threshold_driver, threshold_passenger):
    # profile_dict = profile_to_dict(profile)
    mutation_model = calculate_base_substitution_mutability(profile, cohort_size)
    cohort_size_corrected = cohort_size + 1  # add currently analyzed sample (many samples?) to cohort
    # pprint.pprint(mutation_model)
    results = []
    for mutation_key, mutation_value in mutations_to_rank.items():
        gene, mut = mutation_key
        seq5 = mutation_value['seq5']
        P, Q = mut[0], mut[-1]
        all_substitutions = get_all_codon_substitutions(P, seq5, Q)
        mutability = calculate_codon_mutability(mutation_model, seq5, all_substitutions)

        # print(mutation, mutability)

        # Observed k with pseudocount:
        observed_k = cohort_aa_mutations[gene].get(mut, 0)  # observed in precalculated cohort
        observed_k_pseudocount = observed_k + 1  # we also observe this mutation in analyzed sample

        pval, label, *_ = predict_driver(observed_k_pseudocount, cohort_size_corrected, mutability, threshold_driver, threshold_passenger)
        # print(gene, mut, seq5, all_substitutions, " = ", mutability, observed_k, cohort_size, " = ", pval, label)
        results.append({
            'gene': gene,
            'mutation': mut,
            'mutability': mutability,
            'observed': observed_k_pseudocount,
            'bscore': pval,
            'qvalue': pval,  # temporarily set to pval
            'label': label,
        })
    pvalues = [result['bscore'] for result in results]
    qvalues = []
    if len(pvalues):
        # qvalues = multipletests(pvals=pvalues, method='')[1]
        qvalues = multipletests(pvals=pvalues, method='fdr_bh')[1]
    for i, qvalue in enumerate(qvalues):
        results[i]['qvalue'] = qvalue

    results = list(map(OrderedDict, sorted(results, key=lambda k: k['bscore'])))
    if len(results) == 0:
        return
    df = pd.DataFrame(results, columns=results[0].keys())
    df.drop(df[df.mutability == 0].index, inplace=True)
    try:
        df.to_csv(outfile, sep="\t", index=False, doublequote=False)
    except (BrokenPipeError, IOError):
        pass
