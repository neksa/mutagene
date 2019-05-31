import unittest
import numpy as np
from motifs_in_mutagene import identify_motifs
from mutations_io import read_mutations


class MyTestCase(unittest.TestCase):
    def test_something(self):
        mymotifs = [{
            'name': 'Test N Motif',
            'logo': 'N[T>A]N',
            'motif': 'NTN',
            'position': 1,
            'ref': 'T',
            'alt': 'C',
            "references": "N/A"
        }]
        with open("TCGA-2F-A9KP-01.maf.txt") as f:
            mut = f.read()
            _, mutations_with_context, _ = read_mutations(mut, 'auto', 19)
            observed = identify_motifs(mutations_with_context, mymotifs)
            expected = [{'name': 'Test N Motif', 'motif': 'N[T>A]N', 'enrichment': 0.0, 'pvalue': 1.0, 'mutations': 0}]
            assert observed == expected


if __name__ == '__main__':
    unittest.main()
