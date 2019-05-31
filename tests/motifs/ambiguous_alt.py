import unittest
# import numpy as np
from motifs_in_mutagene import *
from mutations_io import *


class MyTestCase(unittest.TestCase):
    def test_something(self):
        mymotifs = [
            {
                'name': 'Test border Motif',
                'logo': 'T[C>T]W',
                'motif': 'TCW',
                'position': 1,
                'ref': 'C',
                'alt': 'T',
                "references": "N/A"
            },
            {
                'name': 'Test border Motif',
                'logo': 'T[C>T]W',
                'motif': 'TCW',
                'position': 1,
                'ref': 'C',
                'alt': 'G',
                "references": "N/A"
            }
        ]
        with open("TCGA-2F-A9KP-01.maf.txt") as f:
            mut = f.read()
            _, mutations_with_context, _ = read_mutations(mut, 'auto', 19)
            matches = get_enrichment(mutations_with_context, "TCW", 1, "C", "K", 50)
            input_values_1 = get_enrichment(mutations_with_context, mymotifs[0]['motif'], mymotifs[0]['position'],
                                            mymotifs[0]['ref'],
                                            mymotifs[0]['alt'], 50)
            input_values_2 = get_enrichment(mutations_with_context, mymotifs[1]['motif'], mymotifs[1]['position'],
                                            mymotifs[1]['ref'],
                                            mymotifs[1]['alt'], 50)
            expected = (input_values_1[4] + input_values_2[4], input_values_1[5] + input_values_2[5], input_values_1[6],
                        input_values_1[7])
            assert matches[4:] == expected


if __name__ == '__main__':
    unittest.main()
