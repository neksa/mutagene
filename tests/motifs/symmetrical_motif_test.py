import unittest
import numpy as np
from motifs_in_mutagene import identify_motifs, get_enrichment
from mutations_io import read_mutations


class MyTestCase(unittest.TestCase):
    def test_something(self):
        with open("cg_testdata.txt") as f:
            _, mutations_with_context, _ = read_mutations(f.read(), 'auto', 19)
            observed = get_enrichment(mutations_with_context, "CG", 0, "C", "T", 20)
            expected = (0.0, 0, 1.0, 'chi cannot compute due to zero count', 0, 1, 0, 46)
        assert observed == expected


if __name__ == '__main__':
    unittest.main()
