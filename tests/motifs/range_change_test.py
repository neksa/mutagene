import unittest
import numpy as np
from final_maf_code import get_enrichment


class MyTestCase(unittest.TestCase):
    def test_something(self):
        observed = get_enrichment("test_mut.txt", "CTW", 1, "T", "C", 1, assembly=37)  # seq: CCTTG, rev_seq: CAAGG
        expected = (0.0, 0, 0.25, 0.5049850750938457, 1, 0, 0, 3)
        assert observed[4:6] == expected[4:6]  # compare enrichment input


if __name__ == '__main__':
    unittest.main()
