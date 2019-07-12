import unittest
# import numpy as np
from mutagene.motifs import *
# purpose: check that no over counting occurs when overlapping mutations present


class MyTestCase(unittest.TestCase):
    def test_something(self):
        mymotifs = {'motif': 'TCT',
                    'position': 1,
                    'ref': 'C',
                    'alt': 'G'}

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

        observed = get_enrichment(mutations_with_context, mymotifs['motif'], mymotifs['position'], mymotifs['ref'],
                                  mymotifs['alt'], 4, "*")

        assert int(observed['bases_mutated_in_motif']) == 2 \
            and int(observed['bases_not_mutated_in_motif']) == 1 \
            and int(observed['bases_mutated_not_in_motif']) == 0


if __name__ == '__main__':
    unittest.main()
