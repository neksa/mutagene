import unittest
from motifs_in_mutagene import get_enrichment
from mutations_io import read_mutations
import numpy as np


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
            }
        ]
        with open("TCGA-2F-A9KP-01.maf.txt") as f:
            _, mutations_with_context, _ = read_mutations(f.read(), 'auto', 19)
            observed = get_enrichment(mutations_with_context, "TCW", 1, "C", "T", 50)
            expected = (34, 64, 1218, 6927)
            assert observed[4:] == expected


if __name__ == '__main__':
    unittest.main()
