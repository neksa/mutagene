import unittest
from mutagene.motifs import *


class AmbiguousTestCase(unittest.TestCase):
    """purpose: check motif processing of ambiguous bases in alt"""

    def test_ambiguous(self):
        mymotifs = [
            {
                'name': 'Test border Motif',
                'motif': 'AC',
                'position': 1,
                'ref': 'C',
                'alt': 'A',
            },
            {
                'name': 'Test border Motif',
                'motif': 'AC',
                'position': 1,
                'ref': 'C',
                'alt': 'T',
            }
        ]

        mutations_with_context = [
            ('20', 29628279, '+', "C", "A", [
                ('20', 29628277, "C", '+'),
                ('20', 29628278, "A", '+'),
                ('20', 29628279, "C", '+'),
                ('20', 29628280, "G", '+')]),
            ('20', 29628286, '+', "C", "T", [
                ('20', 29628284, "G", '+'),
                ('20', 29628285, "A", '+'),
                ('20', 29628286, "C", '+'),
                ('20', 29628287, "G", '+')])]

        matches, _ = process_mutations(mutations_with_context, "AC", 1, "C", "W", 2, "A")

        input_values_1, _ = process_mutations(
            mutations_with_context, mymotifs[0]['motif'], mymotifs[0]['position'],
            mymotifs[0]['ref'], mymotifs[0]['alt'], 2, "A")

        input_values_2, _ = process_mutations(
            mutations_with_context, mymotifs[1]['motif'], mymotifs[1]['position'],
            mymotifs[1]['ref'], mymotifs[1]['alt'], 2, "A")

        assert matches['bases_mutated_in_motif'] == input_values_1['bases_mutated_in_motif'] + input_values_2['bases_mutated_in_motif']
        assert matches['bases_mutated_not_in_motif'] == input_values_1['bases_mutated_not_in_motif']
        assert matches['bases_mutated_not_in_motif'] == input_values_2['bases_mutated_not_in_motif']
        assert input_values_2['bases_not_mutated_not_in_motif'] == matches['bases_not_mutated_not_in_motif']
        assert input_values_2['bases_not_mutated_not_in_motif'] == input_values_1['bases_not_mutated_not_in_motif']


if __name__ == '__main__':
    unittest.main()
