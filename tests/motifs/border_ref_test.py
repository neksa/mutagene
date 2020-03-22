import unittest
from mutagene.motifs import *


class BorderTestCase(unittest.TestCase):
    """
    purpose: check processing of mutation w/ bases elsewhere in motif
    """

    def test_border_ref(self):
        mymotifs = {
            'motif': 'CTT',
            'position': 1,
            'ref': 'T',
            'alt': 'C'
        }

        mutations_with_context = [
            ('20', 29628279, '+', "T", "C", [
                ('20', 29628278, "C", '+'),
                ('20', 29628279, "T", '+'),
                ('20', 29628280, "T", '+')]),
            ('20', 297, '+', "T", "A", [
                ('20', 296, "C", '+'),
                ('20', 297, "T", '+'),
                ('20', 298, "T", '+')]),
            ('20', 101, '+', "T", "C", [
                ('20', 100, "T", '+'),
                ('20', 101, "T", '+'),
                ('20', 102, "T", '+')]),
            ('20', 11, '-', "A", "G", [
                ('20', 10, "A", '-'),
                ('20', 11, "A", '-'),
                ('20', 12, "G", '-')])]

        observed, saved_matches = process_mutations(
            mutations_with_context,
            mymotifs['motif'],
            mymotifs['position'],
            mymotifs['ref'], mymotifs['alt'], 1, "A")

        assert observed['bases_mutated_in_motif'] == 2
        assert observed['bases_mutated_not_in_motif'] == 1
        assert observed['bases_not_mutated_in_motif'] == 1


if __name__ == '__main__':
    unittest.main()
