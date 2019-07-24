import unittest
from mutagene.motifs import *
# purpose: make sure motifs that only differ by "N' process the same way


class MyTestCase(unittest.TestCase):
    def test_something(self):
        mymotifs = {
            'motif': 'YC',
            'position': 1,
            'ref': 'C',
            'alt': 'T',
        }

        mutations_with_context = [('20', 29628279, '+', "C", "T", [('20', 29628278, "T", '+'),
                                                                    ('20', 29628279, "C", '+'),
                                                                   ('20', 29628280, "A", '+')]),

                                  ('20', 2962, '-', "C", "T", [('20', 2961, "C", '-'),
                                                               ('20', 2962, "C", '-'),
                                                              ('20', 2963, "A", '-')])]

        observed = get_enrichment(mutations_with_context, mymotifs['motif'], mymotifs['position'], mymotifs['ref'], mymotifs['alt'], 1, "*")
        assert int(observed['bases_mutated_in_motif']) == 2


if __name__ == '__main__':
    unittest.main()
