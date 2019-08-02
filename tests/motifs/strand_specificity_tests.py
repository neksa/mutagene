import pytest
import numpy as np
from mutagene.motifs import *


@pytest.mark.parametrize(
    "mutations_with_context",
    [
        [
            ('20', 29628279, '+', 'C', "T",
             [('20', 29628278, 'T', '+'),
              ('20', 29628279, 'C', '+'),
              ('20', 29628280, 'A', '+')])
        ],

    ]
)
def test_strand_for(mutations_with_context):
    # counts for reverse complimentary and forward motifs should be same when testing opposite strands
    observed_TR = process_mutations(mutations_with_context, "YC", 1, "C", "T", 1, "+")
    observed_non_TR = process_mutations(mutations_with_context, "GR", 0, "G", "A", 1, "-")
    assert observed_TR == observed_non_TR

    observed_non_TR_opp = process_mutations(mutations_with_context, "GR", 0, "G", "A", 1, "+")
    observed_TR_opp = process_mutations(mutations_with_context, "YC", 1, "C", "T", 1, "-")
    assert observed_TR_opp == observed_non_TR_opp

    observed = process_mutations(mutations_with_context, "GR", 0, "G", "A", 1, "=")

    assert observed["bases_mutated_in_motif"] == observed_non_TR_opp["bases_mutated_in_motif"] + \
           observed_TR["bases_mutated_in_motif"]

    assert observed["bases_mutated_not_in_motif"] == observed_non_TR_opp["bases_mutated_not_in_motif"] + \
           observed_TR["bases_mutated_not_in_motif"]

    assert observed["bases_not_mutated_in_motif"] == observed_non_TR_opp["bases_not_mutated_in_motif"] + \
           observed_TR["bases_not_mutated_in_motif"]

    assert observed["bases_not_mutated_not_in_motif"] == observed_non_TR_opp["bases_not_mutated_not_in_motif"] + \
           observed_TR["bases_not_mutated_not_in_motif"]


@pytest.mark.parametrize(
    "mutations_with_context",
    [
        [
            ('20', 29628279, '-', 'C', "T",
             [('20', 29628278, 'T', '-'),
              ('20', 29628279, 'C', '-'),
              ('20', 29628280, 'A', '-')])
        ],

    ]
)
def test_strand_rev(mutations_with_context):
    # making sure still works with negative reference strand
    observed_TR = process_mutations(mutations_with_context, "YC", 1, "C", "T", 1, "+")
    observed_non_TR = process_mutations(mutations_with_context, "GR", 0, "G", "A", 1, "-")
    assert observed_TR == observed_non_TR

    observed_non_TR_opp = process_mutations(mutations_with_context, "GR", 0, "G", "A", 1, "+")
    observed_TR_opp = process_mutations(mutations_with_context, "YC", 1, "C", "T", 1, "-")
    assert observed_TR_opp == observed_non_TR_opp

    observed = process_mutations(mutations_with_context, "GR", 0, "G", "A", 1, "=")

    assert observed["bases_mutated_in_motif"] == observed_non_TR_opp["bases_mutated_in_motif"] + \
           observed_TR["bases_mutated_in_motif"]

    assert observed["bases_mutated_not_in_motif"] == observed_non_TR_opp["bases_mutated_not_in_motif"] + \
           observed_TR["bases_mutated_not_in_motif"]

    assert observed["bases_not_mutated_in_motif"] == observed_non_TR_opp["bases_not_mutated_in_motif"] + \
           observed_TR["bases_not_mutated_in_motif"]

    assert observed["bases_not_mutated_not_in_motif"] == observed_non_TR_opp["bases_not_mutated_not_in_motif"] + \
           observed_TR["bases_not_mutated_not_in_motif"]
