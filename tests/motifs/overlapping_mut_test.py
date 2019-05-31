import unittest
import numpy as np
from motifs_in_mutagene import identify_motifs, get_enrichment
from mutations_io import read_mutations


class MyTestCase(unittest.TestCase):
    def test_something(self):
        with open("TCGA-2F-A9KP-01.maf.txt") as f:
            _, mutations, _ = read_mutations(f.read(), 'MAF', 19)
            observed = get_enrichment(mutations, "TCW", 1, "C", "G", 50)
            expected = (23, 10, 1218, 6927)
        assert observed[4:] == expected
        #assert np.isclose([observed[0], expected[0]])  # compare enrichment
        #assert observed[1] == expected[1]  # compare mut burden
        #assert np.isclose([observed[2], expected[2]], [observed[3], expected[3]])  # compare p-values
        #assert observed[4:] == expected[4:]  # compare enrichment input


if __name__ == '__main__':
    unittest.main()
