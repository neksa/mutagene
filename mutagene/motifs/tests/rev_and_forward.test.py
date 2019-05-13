import unittest
from motifs_in_mutagene import identify_motifs, get_enrichment
from mutations_io import read_mutations


class MyTestCase(unittest.TestCase):
    def test_something(self):
        mymotifs = [{
        'name': 'Reverse APOBEC',
        'logo': 'W[G>M]A',
        'motif': 'WGA',
        'position': 1,
        'ref': 'G',
        'alt': 'M',
        'references': 'N/A'
    }]
        with open("TCGA-2F-A9KP-01.maf.txt") as f:
            _, mutations_with_context, _ = read_mutations(f.read(), 'auto', asm=37)
            matches = get_enrichment(mutations_with_context, "TCW", 1,
                        "C", "K", 50)
            rev_match = get_enrichment(mutations_with_context, mymotifs[0]['motif'], mymotifs[0]['position'],
                        mymotifs[0]['ref'], mymotifs[0]['alt'], 50)
            print(matches, rev_match)
            assert matches == rev_match


if __name__ == '__main__':
    unittest.main()
