import pytest
from mutagene.motifs import *


# (record.CHROM, record.POS, record.REF, record.ALT)]
@pytest.mark.parametrize("mutation", [('1', 123, 'A', 'C')])
def test_mutated_base(mutation):
    assert mutated_base(mutation, 'A', 'C')
    assert not mutated_base(mutation, 'A', 'T')


@pytest.mark.parametrize("mutation", [('1', 123, 'A', 'C')])
def test_mutated_base_ambiguous(mutation):
    assert mutated_base(mutation, 'A', 'N')
    assert mutated_base(mutation, 'N', 'C')
    assert mutated_base(mutation, 'N', 'N')
    assert mutated_base(mutation, 'M', 'M')  # M = AC
    assert mutated_base(mutation, 'W', 'S')
    assert mutated_base(mutation, 'R', 'Y')
    assert mutated_base(mutation, 'V', 'H')


def test_scanf():
    custom_motif = 'A[C>T]G'
    s = scanf_motif(custom_motif)
    assert len(s) == 1
    assert s[0]['ref'] == 'C'
    assert s[0]['alt'] == 'T'


@pytest.mark.parametrize(
    "sequence,rev_comp_seq",
    [
        (
            [('chr1', 123, 'A', '+'), ('chr1', 124, 'C', '+'), ('chr1', 125, 'G', '+'), ('chr1', 126, 'T', '+'), ],
            [('chr1', 126, 'A', '-'), ('chr1', 125, 'C', '-'), ('chr1', 124, 'G', '-'), ('chr1', 123, 'T', '-')]
        )
    ])
def test_rev_complementary(sequence, rev_comp_seq):
    for calculated, expected in zip(get_rev_comp_seq(sequence), rev_comp_seq):
        assert calculated == expected


@pytest.mark.parametrize(
    'sequence,motif,motif_position,expected_matches',
    [
        (
            [('chr1', 123, 'A', '+'), ('chr1', 124, 'C', '+'), ('chr1', 125, 'G', '+'), ('chr1', 126, 'T', '+'), ],
            'CG',
            0,
            [('chr1', 124, 'C', '+')]
        ),
        (
            [('chr1', 124, 'C', '+'), ('chr1', 125, 'G', '+')],
            'CG',
            0,
            [('chr1', 124, 'C', '+')]
        ),
        (
            [('chr1', 124, 'C', '-'), ('chr1', 125, 'G', '-')],
            'CG',
            0,
            [('chr1', 124, 'C', '-')]
        ),
        (
            [('chr1', 124, 'C', '+'), ('chr1', 125, 'C', '+')],
            'CG',
            0,
            []
        )
    ],
)
def test_find_matching_motifs(sequence, motif, motif_position, expected_matches):
    matches = list(find_matching_motifs(sequence, motif, motif_position))
    assert matches == expected_matches


@pytest.mark.parametrize(
    'sequence,ref,motif,motif_position,expected_matches',
    [
        (
            [('chr1', 123, 'A', '+'), ('chr1', 124, 'C', '+'), ('chr1', 125, 'G', '+'), ('chr1', 126, 'T', '+'), ],
            'C',
            'CG',
            0,
            [('chr1', 124, 'C', '+')]
        ),
        (
            [('chr1', 124, 'C', '+'), ('chr1', 125, 'G', '+')],
            'C',
            'CG',
            0,
            [('chr1', 124, 'C', '+')]
        ),
        (
            [('chr1', 124, 'C', '-'), ('chr1', 125, 'G', '-')],
            'C',
            'CG',
            0,
            [('chr1', 124, 'C', '-')]
        ),
        (
            [('chr1', 10, 'C', '+'), ('chr1', 11, 'C', '+')],
            'C',
            'CG',
            0,
            [('chr1', 10, 'C', '+')]
        ),
        (
            [('chr1', 10, 'C', '+'), ('chr1', 11, 'C', '+'), ('chr1', 12, 'C', '+'), ('chr1', 13, 'C', '+')],
            'C',
            'CG',
            0,
            [('chr1', 10, 'C', '+'), ('chr1', 11, 'C', '+'), ('chr1', 12, 'C', '+')]
        )
    ],
)
def test_find_matching_bases(sequence, ref, motif, motif_position, expected_matches):
    matches = list(find_matching_bases(sequence, ref, motif, motif_position))
    assert matches == expected_matches
