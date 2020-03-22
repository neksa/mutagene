import unittest
# import numpy as np
from mutagene.motifs import *


class MyTestCase(unittest.TestCase):
    """ purpose: test motif w/ variant base outside mutation """

    def test_something(self):
        mymotifs = {
            'name': 'Test border Motif',
            'logo': 'T[C>T]W',
            'motif': 'TCW',
            'position': 1,
            'ref': 'C',
            'alt': 'T',
            "references": "N/A"
        }

        mutations_with_context = [
            ('20', 101, '+', "C", "T", [
                ('20', 97, "T", '+'),
                ('20', 98, "C", '+'),
                ('20', 99, "T", '+'),

                ('20', 100, "T", '+'),
                ('20', 101, "C", '+'),
                ('20', 102, "T", '+'),

                ('20', 103, "T", '+'),
                ('20', 104, "C", '+'),
                ('20', 105, "T", '+')])]

        mutations_with_context_2 = [
            ('20', 110, '+', "C", "T", [
                ('20', 106, "T", '+'),
                ('20', 107, "C", '+'),
                ('20', 108, "A", '+'),

                ('20', 109, "T", '+'),
                ('20', 110, "C", '+'),
                ('20', 111, "A", '+'),

                ('20', 112, "T", '+'),
                ('20', 113, "C", '+'),
                ('20', 114, "A", '+')])]

        observed_1, saved_matches = process_mutations(
            mutations_with_context,
            mymotifs['motif'],
            mymotifs['position'],
            mymotifs['ref'],
            mymotifs['alt'],
            4,
            "A")

        observed_2, saved_matches = process_mutations(
            mutations_with_context_2,
            mymotifs['motif'],
            mymotifs['position'],
            mymotifs['ref'],
            mymotifs['alt'],
            4,
            "A")

        assert observed_1['bases_mutated_in_motif'] == observed_2['bases_mutated_in_motif']
        assert observed_1['bases_not_mutated_in_motif'] == observed_2['bases_not_mutated_in_motif']
        assert observed_1['bases_mutated_not_in_motif'] == observed_2['bases_mutated_not_in_motif']
        assert observed_1['bases_not_mutated_not_in_motif'] == observed_2['bases_not_mutated_not_in_motif']


if __name__ == '__main__':
    unittest.main()
