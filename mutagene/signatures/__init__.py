
from mutagene.dna import nucleotides


def get_dummy_signatures_lists():
    """
    Generate 6 dummy signatures
    Each will have uniform non-zero frequencies corresponding to one mutation type
    Format them as lists
    """
    dummy_signatures = []
    for mutation in (("C", "A"), ("C", "T"), ("C", "G"), ("T", "A"), ("T", "C"), ("T", "G")):
        values = []
        for p5 in nucleotides:
            for p3 in nucleotides:
                for x in "CT":
                    for y in nucleotides:
                        if x != y:
                            if mutation == (x, y):
                                values.append(1.0 / 16.0)
                            else:
                                values.append(0.0)
        name = mutation[0] + " to " + mutation[1]
        dummy_signatures.append((name, values))
    return dummy_signatures
