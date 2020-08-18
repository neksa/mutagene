
from mutagene.motifs import *
from pprint import pprint

import logging
logger = logging.getLogger(__name__)


def test_enrichment():
    motif_mutation_count = 0
    stat_mutation_count = 8
    stat_motif_count = 0
    stat_ref_count = 5285

    contingency_table = make_contingency_table(
        motif_mutation=motif_mutation_count,
        no_motif_mutation=stat_mutation_count,
        motif_no_mutation=stat_motif_count,
        no_motif_no_mutation=stat_ref_count)

    logger.debug("\n" + contingency_table.to_string() + "\n")

    logger.debug("({} / ({} + {})  ) /  ({} / ({} + {}))".format(
        contingency_table.loc['mutation', 'motif'],
        contingency_table.loc['mutation', 'motif'],
        contingency_table.loc['mutation', 'no motif'],
        contingency_table.loc['no mutation', 'motif'],
        contingency_table.loc['no mutation', 'motif'],
        contingency_table.loc['no mutation', 'no motif']))

    contingency_table = Haldane_correction(contingency_table, pseudocount=1000)

    enrichment = risk_ratio = calculate_RR_for_motif(contingency_table)  # enrichment = risk ratio
    odds_ratio = calculate_OR(contingency_table)

    p_val = get_stats(contingency_table, stat_type=None)

    mut_load = calculate_mutation_load(motif_mutation_count, enrichment)

    result = {
        'enrichment': enrichment,  # AKA risk ratio
        'odds_ratio': odds_ratio,
        'mutation_load': math.ceil(mut_load),
        'pvalue': p_val,
        'bases_mutated_in_motif': motif_mutation_count,
        'bases_mutated_not_in_motif': stat_mutation_count,
        'bases_not_mutated_in_motif': stat_motif_count,
        'bases_not_mutated_not_in_motif': stat_ref_count
    }
    pprint(result)


if __name__ == '__main__':
    test_enrichment()
