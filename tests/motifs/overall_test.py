import unittest
import numpy as np
from motifs_in_mutagene import identify_motifs
from mutations_io import read_mutations


class MyTestCase(unittest.TestCase):
    def test_something(self):
        with open("TCGA-2F-A9KP-01.maf.txt") as f:
            mut = f.read()
            _, mutations_with_context, _ = read_mutations(mut, 'auto', 19)
            observed = identify_motifs(mutations_with_context)
            func = observed[0]
            expected = (4.380675009985355, 44)
        assert np.isclose([func["enrichment"]], [expected[0]])
        assert np.isclose([func["mutations"]], [expected[1]])


if __name__ == '__main__':
    unittest.main()