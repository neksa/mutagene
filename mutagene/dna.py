from itertools import product

nucleotides = "ACGT"  # this order of nucleotides is important for reversing
mutation_contexts = [a + b for a in nucleotides for b in nucleotides]
amino_acids = "ACDEFGHIKLMNPQRSTVWY"
amino_acids_with_stop = amino_acids + "*"
# complementary_nucleotides = reversed(nucleotides)
complementary_nucleotide = dict(zip(nucleotides, reversed(nucleotides)))
complementary_nucleotide['N'] = 'N'

complementary_context = {x: complementary_nucleotide[x[1]] + complementary_nucleotide[x[0]] for x in mutation_contexts}
complementary_trinucleotide = dict(zip(
    ["".join(x) for x in product(nucleotides, repeat=3)],
    ["".join(reversed(x)) for x in product(reversed(nucleotides), repeat=3)]))

bases_dict = {
    "A": "A", "G": "G", "T": "T", "C": "C",
    "W": "AT", "S": "CG", "M": "AC", "K": "GT", "R": "AG", "Y": "CT",
    "B": "TCG", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ATGC"}

extended_nucleotides = "ACTGWSMKRYBDHVN"
complementary_extended_nucleotide = dict(zip(extended_nucleotides, "TGACWSKMYRVHDBN"))

comp_dict = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "W": "AT", "S": "CG", "K": "AC", "M": "GT", "Y": "AG", "R": "CT",
    "V": "TCG", "H": "AGT", "D": "ACT", "B": "ACG", "N": "ATGC"}

codon_table = {
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "AAA": "K", "AAG": "K",
    "AAT": "N", "AAC": "N",
    "ATG": "M",
    "GAT": "D", "GAC": "D",
    "TTT": "F", "TTC": "F",
    "TGT": "C", "TGC": "C",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "GAA": "E", "GAG": "E",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "TGG": "W",
    "CAT": "H", "CAC": "H",
    "TAT": "Y", "TAC": "Y",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TAG": "*", "TGA": "*", "TAA": "*"
}

exome_trinucleotides = {
    "GCA": 1870205, "ACT": 1317607, "GCC": 2018826, "CCT": 2045943, "GTC": 1272257, "ATT": 1226179, "CTC": 1989259,
    "ACA": 1614930, "ATC": 1343451, "ACG": 602918, "TTC": 1901940, "GTT": 1149725, "GCG": 857006, "GTG": 1685105,
    "ACC": 1369960, "CCA": 2359526, "TTG": 1588902, "ATA": 841828, "TCA": 1853413, "CCG": 1009679, "TTA": 774505,
    "TCG": 640487, "ATG": 1654761, "GTA": 798348, "CTT": 1881403, "GCT": 1983552, "CTA": 713831, "TTT": 1756413,
    "CCC": 1827705, "TCC": 2026380, "TCT": 2000322, "CTG": 2769315}

aa_short = amino_acids_with_stop
aa_long = ["Ala", "Leu", "Pro", "Gly", "Met", "Ser", "Thr", "Trp", "Ile", "Val", "Cys", "Asp", "Glu", "Phe", "His", "Lys", "Asn", "Gln", "Arg", "Tyr", "STOP"]
aa_dict = dict(zip(aa_short, aa_long))

bases_dict = {"A": "A", "G": "G", "T": "T", "C": "C", "W": "AT", "S": "CG", "M": "AC", "K": "GT", "R": "AG", "Y": "CT",
              "B": "TCG", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ATGC"}

comp_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "W": "AT", "S": "CG", "K": "AC", "M": "GT", "Y": "AG",
             "R": "CT", "V": "TCG", "H": "AGT", "D": "ACT", "B": "ACG", "N": "ATGC"}

chromosome_name_mapping = {
    "chr23": "chrX",
    "chr24": "chrY",
    "chr25": "chrXY",
    "chr26": "chrM",
}
