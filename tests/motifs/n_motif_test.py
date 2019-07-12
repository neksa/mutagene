import unittest
from mutagene.motifs import *
# purpose: check  processing of motifs w/ "N"


class MyTestCase(unittest.TestCase):
    def test_something(self):
        mymotifs = {
            'motif': 'NTN',
            'position': 1,
            'ref': 'T',
            'alt': 'C',
        }

        mutations_with_context = [('20', 29628279, '+', "T", "C", [('20', 29628278, "C", '+'),
                                                                    ('20', 29628279, "T", '+'),
                                                                    ('20', 29628280, "T", '+')]),

                                  ('20', 11, '-', "A", "G", [('20', 10, "A", '-'),
                                                               ('20', 11, "A", '-'),
                                                               ('20', 12, "G", '-')])
                                  ]

        observed = get_enrichment(mutations_with_context, mymotifs['motif'], mymotifs['position'], mymotifs['ref'], mymotifs['alt'], 1, "*")
        assert int(observed['bases_mutated_in_motif']) == 2


if __name__ == '__main__':
    unittest.main()
