import pandas as pd

from collections import defaultdict, OrderedDict
from functools import lru_cache, reduce

from scipy.stats import binom_test
from statsmodels.stats.multitest import multipletests

from .dna import *
from .io import get_signature_attributes_dict

import logging
logger = logging.getLogger(__name__)


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
def predict_driver(observed, N, p):
    """
    Binomial P-value based on the number of observed mutations, number of samples and expected mutability
    """
    # TODO: thresholds should be defined in the config
    # need more tumor samples than mutations for sure

    if N > 0 and observed > 0 and observed <= N:
        pvalue = binom_test(observed, N, p * 1e-6, alternative='greater')
    else:
        pvalue = 1.0

    if pvalue <= 8.030647e-05:  # ? Max MCC
        prediction = "Driver"
        prediction_value = 1
    elif pvalue <= 0.003440945:  # ? 10 % FPR
        prediction = "Potential driver"
        prediction_value = 0.5
    elif pvalue < 1.0:
        prediction = "Passenger"
        prediction_value = 0
    else:
        prediction = "Undefined"
        prediction_value = -1

    return pvalue, prediction, prediction_value


def profile_to_dict(counts_profile):
    mutation_prob = defaultdict(dict)
    total_count = sum(counts_profile)
    freqs = [x / total_count for x in counts_profile]

    for i, attrib in enumerate(get_signature_attributes_dict()):
        x, y = attrib['mutation']
        p5, p3 = attrib['context']
        tri = p5 + x + p3
        mutation_prob[tri][y] = freqs[i]
        mutation_prob[complementary_trinucleotide[tri]][complementary_nucleotide[y]] = freqs[i]
    return mutation_prob


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


def rank(mutations_to_rank, outfile, profile, cohort_aa_mutations, cohort_size):
    mutation_model = profile_to_dict(profile)
    results = []
    for gene, mut, seq5 in mutations_to_rank:
        P, Q = mut[0], mut[-1]
        all_substitutions = get_all_codon_substitutions(P, seq5, Q)
        mutability = calculate_codon_mutability(mutation_model, seq5, all_substitutions)

        # Observed N with pseudocount:
        observed_k = cohort_aa_mutations[gene].get(mut, 0) + 1
        N = cohort_size + 1

        pval, label, *_ = predict_driver(observed_k, N, mutability)
        # print(gene, mut, seq5, all_substitutions, " = ", mutability, observed_k, cohort_size, " = ", pval, label)
        results.append({
            'gene': gene,
            'mutation': mut,
            'mutability': mutability,
            'bscore': pval,
            'qvalue': pval,  # temporarily set to pval
            # 'label': label,
        })
    pvalues = [result['bscore'] for result in results]
    qvalues = []
    if len(pvalues):
        qvalues = multipletests(pvals=pvalues, method='fdr_bh')[1]
    for i, qvalue in enumerate(qvalues):
        results[i]['qvalue'] = qvalue

    results = list(map(OrderedDict, sorted(results, key=lambda k: k['qvalue'])))
    if len(results) == 0:
        return
    df = pd.DataFrame(results, columns=results[0].keys())
    df.drop(df[df.mutability == 0].index, inplace=True)
    df.to_csv(outfile, sep="\t", index=False, doublequote=False)
