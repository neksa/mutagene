import pytest
from mutagene.motifs import *


@pytest.mark.parametrize(
    "motif, position, alt, window_size, strand, mutations_with_context",
    [(
        'YC',
        1,
        'T',
        1,
        'A',
        [
            ('20', 29628279, '+', 'C', "T", [
                ('20', 29628278, 'T', '+'),
                ('20', 29628279, 'C', '+'),
                ('20', 29628280, 'A', '+')]),

            ('20', 2962, '-', 'C', "T", [
                ('20', 2961, 'C', '-'),
                ('20', 2962, 'C', '-'),
                ('20', 2963, 'A', '-')])])]
)
def test_asymmetrical_rev_motif(motif, position, alt, window_size, strand, mutations_with_context):
    # purpose: make sure motifs that only differ by "N' process the same way
    observed, _ = process_mutations(mutations_with_context, motif, position, motif[position], alt, window_size, strand)
    assert observed['bases_mutated_in_motif'] == 2
