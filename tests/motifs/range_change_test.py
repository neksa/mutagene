import unittest
# import numpy as np
from mutagene.motifs import *
# purpose: check processing of diff range sizes


class MyTestCase(unittest.TestCase):
    def test_something(self):
        mymotifs = { 'motif': 'TCT',
            'position': 1,
            'ref': 'C',
            'alt': 'G'}

        mutations_with_long_range = [('20', 101, '+', "C", "G",
                                                             [('20', 97, "T", '+'),
                                                              ('20', 98, "C", '+'),
                                                              ('20', 99, "T", '+'),

                                                              ('20', 100, "T", '+'),
                                                              ('20', 101, "C", '+'),
                                                              ('20', 102, "T", '+'),

                                                              ('20', 103, "T", '+'),
                                                              ('20', 104, "C", '+'),
                                                              ('20', 105, "T", '+')
                                                              ])]

        mutations_with_short_range = [('20', 101, '+', "C", "G",
                                                              [('20', 100, "T", '+'),
                                                              ('20', 101, "C", '+'),
                                                              ('20', 102, "T", '+')
                                                              ])]


        observed_long = get_enrichment(mutations_with_long_range, mymotifs['motif'], mymotifs['position'],
                                       mymotifs['ref'], mymotifs['alt'], 4, "=")

        observed_short = get_enrichment(mutations_with_short_range, mymotifs['motif'], mymotifs['position'],
                                       mymotifs['ref'], mymotifs['alt'], 1, "=")

        print(observed_long)
        print(observed_short)

        assert int(observed_long['bases_mutated_in_motif']) == int(observed_short['bases_mutated_in_motif']) \
            and int(observed_long['bases_mutated_not_in_motif']) == int(observed_short['bases_mutated_not_in_motif']) \
            and int(observed_long['bases_not_mutated_in_motif']) == int(observed_short['bases_not_mutated_in_motif'] + 2)


if __name__ == '__main__':
    unittest.main()
