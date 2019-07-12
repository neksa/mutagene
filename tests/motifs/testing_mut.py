import unittest
import numpy as np
from final_maf_code import get_enrichment


def testing_mut(self):
    observed = get_enrichment("CGA-2F-A9KP-01.maf.txt", "TCW", 1, "C", "K", 50, assembly=37)
    expected = (4.394267654751525,
                45,
                8.041476308812957e-15,
                6.379300415739294e-19,
                57,
                74,
                1240,
                7074)
    assert np.isclose([observed[0], expected[0]]) #compare enrichment
    assert observed[1] == expected[1] #compare mut burden
    assert np.isclose([observed[2], expected[2]], [observed[3], expected[3]]) #compare p-values
    assert observed[4:] == expected[4:] #compare enrichment input

    observed1 = get_enrichment("test_mutations.txt", "TCW", 1, "C", "K", 50, assembly=37)
    expected1 = (2.284126365054602,
                 10,
                 0.005675513210048566,
                 0.006578784725027336,
                 17,
                 40,
                 641,
                 3445)
    assert np.isclose([observed1[0], expected1[0]]) #compare enrichment
    assert observed1[1] == expected1[1] #compare mut burden
    assert np.isclose([observed1[2], expected1[2]], [observed1[3], expected1[3]]) #compare p-values
    assert observed1[4:] == expected1[4:] #compare enrichment input

    observed2 = get_enrichment("test_mut.txt", "CTT", 1, "T", "C", 1, assembly=37)
    expected2 = (0.0, 0, 0.2500000000000001, 0.5049850750938457, 1, 0, 0, 3)
    assert np.isclose([observed2[0], expected2[0]])  # compare enrichment
    assert observed2[1] == expected2[1]  # compare mut burden
    assert np.isclose([observed2[2], expected2[2]], [observed2[3], expected2[3]])  # compare p-values
    assert observed2[4:] == expected2[4:]  # compare enrichment input

    observed3 = get_enrichment("test_mut.txt", "TCA", 1, "C", "A", 1, assembly=37)
    expected3 = (0.0, 0, 0.5000000000000003, 1.0, 1, 0, 0, 1)
    assert np.isclose([observed3[0], expected3[0]])  # compare enrichment
    assert observed3[1] == expected3[1]  # compare mut burden
    assert np.isclose([observed3[2], expected3[2]], [observed3[3], expected3[3]])  # compare p-values
    assert observed3[4:] == expected3[4:]  # compare enrichment input


if __name__ == '__main__':

    unittest.main()
