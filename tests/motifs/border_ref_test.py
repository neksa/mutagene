import unittest
import numpy as np
from motifs_in_mutagene import identify_motifs
from mutations_io import read_mutations


class MyTestCase(unittest.TestCase):
    def test_something(self):
        mymotifs = [{
            'name': 'Test border Motif',
            'logo': 'C[T>C]T',
            'motif': 'CTT',
            'position': 1,
            'ref': 'T',
            'alt': 'C',
            "references": "N/A"
        }]

        with open("test_mut_pycharm.txt") as f:
            mut = f.read()
            _, mutations_with_context, _ = read_mutations(mut, 'auto', 19)
            observed = identify_motifs(mutations_with_context, mymotifs)
            expected = (0.0, 0, 0.08547008547008556, 'chi cannot compute due to zero count', 1, 0, 9, 107)
            print(observed)
            expected = ""
            assert observed == expected


if __name__ == '__main__':
    unittest.main()
