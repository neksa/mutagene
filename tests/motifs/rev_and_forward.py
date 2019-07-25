import unittest
# import numpy as np
from mutagene.motifs import *
# purpose: check processing of reverse and forward motifs


class MyTestCase(unittest.TestCase):
    def test_something(self):
        mymotifs_for = {
        'motif': 'TCT',
        'position': 1,
        'ref': 'C',
        'alt': 'T',
    }

        mymotifs_rev = {
            'motif': 'AGA',
            'position': 1,
            'ref': 'G',
            'alt': 'A',
        }

        mutations_with_context = [
                                  ('20', 101, '+', "C", "G", [('20', 97, "T", '+'),
                                                              ('20', 98, "C", '+'),
                                                              ('20', 99, "T", '+'),

                                                              ('20', 100, "T", '+'),
                                                              ('20', 101, "C", '+'),
                                                              ('20', 102, "T", '+'),

                                                              ('20', 103, "T", '+'),
                                                              ('20', 104, "C", '+'),
                                                              ('20', 105, "T", '+')
                                                              ]),

                                  ('20', 104, '+', "C", "G", [('20', 100, "T", '+'),
                                                              ('20', 101, "C", '+'),
                                                              ('20', 102, "T", '+'),

                                                              ('20', 103, "T", '+'),
                                                              ('20', 104, "C", '+'),
                                                              ('20', 105, "T", '+'),

                                                              ('20', 106, "A", '+'),
                                                              ('20', 107, "A", '+'),
                                                              ('20', 108, "A", '+'),
                                                              ])
                                  ]

        observed_for = get_enrichment(mutations_with_context, mymotifs_for['motif'], mymotifs_for['position'],
                                      mymotifs_for['ref'], mymotifs_for['alt'], 4, "=")

        observed_rev = get_enrichment(mutations_with_context, mymotifs_rev['motif'], mymotifs_rev['position'],
                                      mymotifs_rev['ref'], mymotifs_rev['alt'], 4, "=")
        assert int(observed_for['bases_mutated_in_motif']) == int(observed_rev['bases_mutated_in_motif']) \
            and int(observed_for['bases_not_mutated_in_motif']) == int(observed_rev['bases_not_mutated_in_motif']) \
            and int(observed_for['bases_mutated_not_in_motif']) == int(observed_rev['bases_mutated_not_in_motif']) \
            and int(observed_for['bases_not_mutated_not_in_motif']) == int(observed_rev['bases_not_mutated_not_in_motif'])


if __name__ == '__main__':
    unittest.main()
