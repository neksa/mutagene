import logging

from mutagene.motifs import *

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
        no_motif_no_mutation=stat_ref_count,
    )

    logger.debug("\n" + contingency_table.to_string() + "\n")

    logger.debug(
        "({} / ({} + {})  ) /  ({} / ({} + {}))".format(
            contingency_table.loc["mutation", "motif"],
            contingency_table.loc["mutation", "motif"],
            contingency_table.loc["mutation", "no motif"],
            contingency_table.loc["no mutation", "motif"],
            contingency_table.loc["no mutation", "motif"],
            contingency_table.loc["no mutation", "no motif"],
        )
    )

    contingency_table = Haldane_correction(contingency_table, pseudocount=1000)

    enrichment = calculate_RR_for_motif(contingency_table)  # enrichment = risk ratio
    calculate_OR(contingency_table)

    get_stats(contingency_table, stat_type=None)

    calculate_mutation_load(motif_mutation_count, enrichment)


if __name__ == "__main__":
    test_enrichment()
