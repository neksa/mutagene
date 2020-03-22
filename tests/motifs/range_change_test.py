import unittest
from mutagene.motifs import *


class MyTestCase(unittest.TestCase):
    """purpose: check processing of diff range sizes"""

    def test_something(self):
        mymotifs = {
            'motif': 'TCT',
            'position': 1,
            'ref': 'C',
            'alt': 'G'}

        mutations_with_long_range = [
            ('20', 101, '+', "C", "G", [
                ('20', 97, "T", '+'),
                ('20', 98, "C", '+'),
                ('20', 99, "T", '+'),

                ('20', 100, "T", '+'),
                ('20', 101, "C", '+'),
                ('20', 102, "T", '+'),

                ('20', 103, "T", '+'),
                ('20', 104, "C", '+'),
                ('20', 105, "T", '+')])]

        mutations_with_short_range = [
            ('20', 101, '+', "C", "G", [
                ('20', 100, "T", '+'),
                ('20', 101, "C", '+'),
                ('20', 102, "T", '+')])]

        observed_long, saved_matches = process_mutations(
            mutations_with_long_range,
            mymotifs['motif'],
            mymotifs['position'],
            mymotifs['ref'],
            mymotifs['alt'],
            4,
            "A")

        observed_short, saved_matches = process_mutations(
            mutations_with_short_range,
            mymotifs['motif'],
            mymotifs['position'],
            mymotifs['ref'],
            mymotifs['alt'],
            1,
            "A")

        assert observed_long['bases_mutated_in_motif'] == observed_short['bases_mutated_in_motif']
        assert observed_long['bases_mutated_not_in_motif'] == observed_short['bases_mutated_not_in_motif']
        assert observed_long['bases_not_mutated_in_motif'] == observed_short['bases_not_mutated_in_motif'] + 2


if __name__ == '__main__':
    unittest.main()
