import twobitreader as tbr
from collections import defaultdict
from pprint import pprint
from statsmodels.stats.multitest import multipletests

import scipy.stats as stats
import csv
import numpy as np
from itertools import cycle
import matplotlib.pyplot as plt
import math

nucleotides = "ACGT"  # this order of nucleotides is important for reversing
complementary_nucleotide = dict(zip(nucleotides, reversed(nucleotides)))
complementary_nucleotide['N'] = 'N'

TWOBIT_GENOMES_PATH = '/net/pan1/mutagene/data/genomes/'


bases_dict = {"A": "A", "G": "G", "T": "T", "C": "C", "W": "AT", "S": "CG", "M": "AC", "K": "GT", "R": "AG", "Y": "CT",
              "B": "TCG", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ATGC"}

comp_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "W": "AT", "S": "CG", "K": "AC", "M": "GT", "Y": "AG",
             "R": "CT", "V": "TCG", "H": "AGT", "D": "ACT", "B": "ACG", "N": "ATGC"}


def get_stats(motif_mutation_count, mutation_count, motif_count, ref_count):
    p_val = stats.fisher_exact([[mutation_count, motif_mutation_count], [ref_count, motif_count]], alternative = "less")
    p_value = p_val[1]
    chi_array = np.array([[mutation_count, motif_mutation_count], [ref_count, motif_count]])
    chi = stats.chi2_contingency(chi_array)[1]
    # if p_value <= 0.05:
    #  qvalues = multipletests(pvals=p_value, method='fdr_bh')
    #  print(qvalues)
    #  if qvalues[3] <= 0.05:
    #     print("significant")
    # print("odds_ratio: ", "p-value")
    return p_val[0], p_value, chi


def motif_contexts(mutation, motif, motif_position, range_size, assembly=None):
    """
    :param mutation: [(record.CHROM, record.POS, record.REF, record.ALT)]
    :param motif: contains the ref nucleotide and surrounding nucleotides
    :param range_size: number of nucleotides the sequence extends to excluding the motif on each side
    :param assembly: specific genome sequence, default 38
    :return: context of a mutation, range_size included
    """

    genomes_path = TWOBIT_GENOMES_PATH
    twobit_files = {
        38: 'hg38',
        37: 'hg19',
        19: 'hg19'
    }

    chromosome_name_mapping = {
        "chr23": "chrX",
       "chr24": "chrY",
        "chr25": "chrXY",
        "chr26": "chrM",
    }

    if assembly is None:
        assembly = 38

    assert assembly in twobit_files

    twobit_file = genomes_path + "/" + twobit_files[assembly] + ".2bit"
    f = tbr.TwoBitFile(twobit_file)
    cn = complementary_nucleotide
    chrom = str(mutation[0])
    start = mutation[1] - 1  # zero-based numbering; start=0
    chromosome = chrom if chrom.startswith('chr') else 'chr' + chrom
    chromosome = chromosome_name_mapping.get(chromosome, chromosome)
    if chromosome in f:
        try:
            seq = f[chromosome][start - range_size - motif_position:
                                start + len(motif) + range_size - motif_position]
            # print(motif_position)
            seq = seq.upper()
        except Exception as e:
            print("TwoBit exception", str(e), mutation)
    else:
        print("NO CHROM", chromosome)

    return list(zip(
        cycle([mutation[0]]),
        range(mutation[1] - range_size - motif_position, mutation[1] + len(motif) + range_size - motif_position),
        seq,
        cycle("+")))


def read_MAF(filename):
    """
    :param filename: name of MAF file in quotes
    :param motif: contains the ref nucleotide and surrounding nucleotides
    :param motif_position: zero-based numbering, which nucleotide is being mutated in motif
    :param ref: the nucleotide base pre-mutation
    :param alt: the nucleotide base post-mutation
    :param range_size: the number of nucleotides the sequence extends to excluding the motif on each side
    :param assembly: specific genome sequence, default 38
    :return: numbers needed to calc enrichment
    """
    # function counts number of mutations w/o regard for motif (identical to VCF list_murations())
    # cn = complementary_nucleotide
    # mutations = defaultdict(float)
    # N_skipped = 0

    sample_list = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            if line.startswith("Hugo_Symbol"):
                continue
            if len(line) < 10:
                continue

            col_list = line.split("\t")
            if len(col_list) < 13:
                continue

            try:
                # assembly_build = col_list[3]  # MAF ASSEMBLY
                # strand = col_list[7]    # MAF STRAND
                # ID = col_list[2]

                # chromosome is expected to be one or two number or one letter
                chrom = col_list[4]  # MAF CHROM
                if chrom.lower().startswith("chr"):
                    chrom = chrom[3:]
                # if len(chrom) == 2 and chrom[1] not in "0123456789":
                #     chrom = chrom[0]

                pos = int(col_list[5])  # MAF POS START
                pos_end = int(col_list[6])  # MAF POS END
                x = col_list[10]  # MAF REF
                #x = col_list[10]
                #y1 = col_list[11]  # MAF ALT1
                y = col_list[12]
                #y2 = col_list[12]  # MAF ALT2
                sample = col_list[15]  # MAF VARIANT_TYPE
                #sample = col_list[15]

            except:
                # raise
                continue

            if pos != pos_end:
                continue

            # skip if found unexpected nucleotide characters
            #if len(set([x, y1, y2]) - set(nucleotides)) > 0:
                #continue

            if len(set([x, y]) - set(nucleotides)) > 0:
                continue
            # y = y1 if y1 != x else None
            # y = y2 if y2 != x else y

            if y is None:
                continue

            sample_list.append((chrom, pos, x, y))

        return sample_list


def get_rev_comp_seq(sequence):
    # rev_comp_seq = "".join([complementary_nucleotide[i] for i in reversed(sequence)])
    rev_comp_seq = [(i[0], i[1], complementary_nucleotide[i[2]], "-") for i in reversed(sequence)]
    return rev_comp_seq


def mutated_base(mutation, ref, alt):
    """
    :param mutation: [(record.CHROM, record.POS, record.REF, record.ALT)]
    :param ref: the nucleotide base pre-mutation
    :param alt: the nucleotide base post-mutation
    :return: how many mutations match specified ref and alt
    """
    assert ref != alt, "mutation should have different ref and alt nucleotides"

    # makes sure single base substitution
    if mutation[3] and mutation[2] and len(mutation[2]) == 1 \
            and len(mutation[3]) == 1 and len(mutation[3]) == 1 \
            and len(mutation[2]) == 1 and mutation[2] != mutation[3]:
        # mutation matches the substitution
        if mutation[2] in ref and mutation[3] in alt:
            return True


def find_matching_motifs(seq, motif, motif_position):
    # print("Looking for motif {} in {}, {}".format(motif, sequence, len(sequence) - len(motif)))
    for i in range(len(seq) - len(motif) + 1):
        s = seq[i: i + len(motif)]
        for j, char in enumerate(motif):
            if s[j][2] not in bases_dict[char]:
                break
        else:
            yield seq[i+motif_position]


def find_matching_bases(seq, ref, motif, motif_position):
    # print("Looking for motif {} in {}, {}".format(motif, sequence, len(sequence) - len(motif)))
    for i in range(len(seq)):
        s = seq[i]
        if s[2] in ref:
            yield seq[i]


def get_enrichment(filename, motif, motif_position, ref, alt, range_size, assembly=None):
    mutations = read_MAF(filename)
    assert range_size >= 0
    assert len(ref) == 1
    assert len(alt) == 1
    assert 0 <= motif_position < len(motif)
    assert motif[motif_position] == ref

    matching_bases = set()
    matching_motifs = set()
    matching_mutated_motifs = set()
    matching_mutated_bases = set()

    matching_mutated_ref = set()

    # extra loop for sample in sample list
    for mutation in mutations:
        # extract the longest sequence we would ever need (motif + range_size)

        seq = motif_contexts(mutation, motif, motif_position, range_size, assembly)
        rev_seq = get_rev_comp_seq(seq)
        rev_motif_pos = len(motif) - motif_position - 1

        if mutation[2] == bases_dict[ref] or mutation[2] == comp_dict[ref]:
            matching_mutated_ref.add(mutation[0:2])

        # not mutated:
        ref_seq = seq[motif_position:len(seq) - (len(motif) - motif_position) + 1]

        for ref_match in find_matching_bases(ref_seq, ref, motif, motif_position):
            matching_bases.add(ref_match[0:2])

        for motif_match in find_matching_motifs(seq, motif, motif_position):
            matching_motifs.add(motif_match[0:2])

        # rev compl: not mutated:
        ref_rev_seq = rev_seq[rev_motif_pos: len(seq) - (len(motif) - rev_motif_pos) + 1]

        for ref_match in find_matching_bases(ref_rev_seq, ref, motif, len(motif) - motif_position - 1):
            matching_bases.add(ref_match[0:2])

        for motif_match in find_matching_motifs(rev_seq, motif, rev_motif_pos):
            matching_motifs.add(motif_match[0:2])

        # mutated:
        if mutated_base(mutation, bases_dict[ref], bases_dict[alt]):
            matching_mutated_bases.add(mutation[0:2])

            context_of_mutation = seq[range_size: range_size + len(motif)]
            for motif_match in find_matching_motifs(context_of_mutation, motif, motif_position):
                matching_mutated_motifs.add(motif_match[0:2])

        if mutated_base(mutation, comp_dict[ref], comp_dict[alt]):
            matching_mutated_bases.add(mutation[0:2])

            # rev comp:
            context_of_mutation = rev_seq[range_size: range_size + len(motif)]
            for motif_match in find_matching_motifs(context_of_mutation, motif, rev_motif_pos):
                matching_mutated_motifs.add(motif_match[0:2])

    motif_mutation_count = len(matching_mutated_motifs)
    mutation_count = len(matching_mutated_bases - matching_mutated_motifs)

    motif_count = len(matching_motifs - matching_mutated_ref)
    ref_count = len(matching_bases - matching_motifs - matching_mutated_ref)

    stat_motif_count = len(matching_motifs - matching_mutated_motifs - matching_mutated_ref)
    stat_ref_count = len(matching_bases - matching_motifs - matching_mutated_bases - matching_mutated_ref)

    try:
        enrichment = (motif_mutation_count / mutation_count) / (motif_count / ref_count)
    except ZeroDivisionError:
        enrichment = 0.0

    p_val = get_stats(motif_mutation_count, mutation_count, stat_motif_count, stat_ref_count)

    if enrichment > 1 and p_val[1] < 0.05 and p_val[2] < 0.05:
        mut_load = (motif_mutation_count * (enrichment - 1)) / enrichment
    else:
        mut_load = 0.0
    return enrichment, math.ceil(mut_load), p_val[1], p_val[2], motif_mutation_count, mutation_count, motif_count, \
           ref_count


# def get_values(filename, motif, motif_position, ref, alt, assembly=None):
#     mut_burden_count = []
#     range_sizes = []
#     with open("APOBEC_window.csv", "w") as csvfile: #create csv file
#         filewriter = csv.writer(csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL)
#         filewriter.writerow([motif])
#         filewriter.writerow(["range_size", "enrichment", "motif_mut_count", "mut_count", "motif_count", "ref_count", "mutation load", "fisher p-value", "chi p-value"])
#         for val in range(0, 200 + 1, 5):
#             range_sizes.append(val * 2 + len(motif))
#             analysis = get_enrichment(filename, motif, motif_position, ref, alt, val, assembly=assembly)
#             result = analysis[0]
#             mut_burden = analysis[1]
#             mut_burden_count.append(mut_burden)
#             pvalue = analysis[2]
#             chi_pvalue = analysis[3]
#             motif_mut = analysis[4]
#             mut = analysis[5]
#             motif_count = analysis[6]
#             ref_count = analysis[7]
#             filewriter.writerow([val, result, motif_mut, mut, motif_count, ref_count, mut_burden, pvalue, chi_pvalue])
#         plt.plot(range_sizes, mut_burden_count, color = 'b')
#         plt.xlabel("Window Size")
#         plt.ylabel("Mutational Burden")
#         plt.title("APOBEC Mutational Burden Depends Upon Window Size")
#         plt.show()
#         return list(zip(range_sizes, mut_burden_count))


def get_values(filename, motif, motif_position, ref, alt, assembly=None):
    mut_burden_count = []
    range_sizes = []
    with open("APOBEC5_window.csv", "w") as csvfile: #create csv file
        filewriter = csv.writer(csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow([motif])
        filewriter.writerow(["range_size", "enrichment", "motif_mut_count", "mut_count", "motif_count", "ref_count", "mutation load", "fisher p-value", "chi p-value"])
        for val in range(0, 2000 + 1, 100):
            range_sizes.append(val * 2 + len(motif))
            analysis = get_enrichment(filename, motif, motif_position, ref, alt, val, assembly=assembly)
            result = analysis[0]
            mut_burden = analysis[1]
            mut_burden_count.append(mut_burden)
            pvalue = analysis[2]
            chi_pvalue = analysis[3]
            motif_mut = analysis[4]
            mut = analysis[5]
            motif_count = analysis[6]
            ref_count = analysis[7]
            filewriter.writerow([val, result, motif_mut, mut, motif_count, ref_count, mut_burden, pvalue, chi_pvalue])
        plt.plot(range_sizes, mut_burden_count, color='b')
        plt.xlabel("Window Size")
        plt.ylabel("Mutational Burden")
        plt.title("APOBEC Mutational Burden Depends Upon Window Size")
        plt.show()
        return list(zip(range_sizes, mut_burden_count))


if __name__ == '__main__':
    # pprint(get_enrichment("tcga_131.txt", "TCW", 1, "C", "K", 50, assembly=37))
    # pprint(get_enrichment("TCGA-2F-A9KP-01.maf.txt", "WGA", 1, "G", "M", 19, assembly=37))

    pprint(get_enrichment("TCGA-2F-A9KQ-01.maf.txt", "TCW", 1, "C", "K", 50, assembly=37))
    pprint(get_enrichment("TCGA-2F-A9KQ-01.maf.txt", "WGA", 1, "G", "M", 50, assembly=37))

    # pprint(get_enrichment("TCGA-2F-A9KO-01.maf", "YC", 1, "C", "T", 50, assembly=37))

    # pprint(get_enrichment("TCGA-2F-A9KO-01.maf", "NTW", 1, "T", "G", 50, assembly=37))
    # pprint(get_enrichment("TCGA-2F-A9KO-01.maf", "NTW", 1, "T", "C", 50, assembly=37))

    # pprint(get_enrichment("TCGA-2F-A9KO-01.maf", "CTG", 1, "T", "A", 50, assembly=37))

    # pprint(get_enrichment("TCGA-2F-A9KO-01.maf", "CT", 0, "C", "A", 50, assembly=37))
    # pprint(get_enrichment("TCGA-2F-A9KO-01.maf", "CG", 0, "C", "T", 50, assembly=37))

    # pprint(get_enrichment("TCGA-2F-A9KO-01.maf", "NTW", 1, "T", "G", 50, assembly=37))
    # pprint(get_enrichment("TCGA-2F-A9KO-01.maf", "NTW", 1, "T", "C", 50, assembly=37))

    # pprint(get_enrichment("TCGA-2F-A9KO-01.maf", "CG", 0, "C", "A", 50, assembly=37))
    # pprint(get_enrichment("TCGA-2F-A9KO-01.maf", "CG", 0, "C", "T", 50, assembly=37))
    # pprint(get_enrichment("TCGA-2F-A9KO-01.maf", "CG", 0, "C", "G", 50, assembly=37))

