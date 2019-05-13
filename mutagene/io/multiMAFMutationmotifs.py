import twobitreader as tbr
from collections import defaultdict
from pprint import pprint
from statsmodels.stats.multitest import multipletests
import scipy.stats as stats
import csv
import numpy as np
from itertools import cycle
import matplotlib.style as style
import matplotlib.pyplot as plt

#like MAFMutationmotifs but works for many samples

nucleotides = "ACGT"  # this order of nucleotides is important for reversing
complementary_nucleotide = dict(zip(nucleotides, reversed(nucleotides)))
complementary_nucleotide['N'] = 'N'

TWOBIT_GENOMES_PATH = '/net/pan1/mutagene/data/genomes/'


bases_dict = {"A": "A", "G": "G", "T": "T", "C": "C", "W": "AT", "S": "CG", "M": "AC", "K": "GT", "R": "AG", "Y": "CT",
              "B": "TCG", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ATGC"}

comp_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "W": "AT", "S": "CG", "K": "AC", "M": "GT", "Y": "AG",
             "R": "CT",
             "V": "TCG", "H": "AGT", "D": "ACT", "B": "ACG", "N": "ATGC"}


def get_stats(motif_mutation_count, mutation_count, motif_count, ref_count):
    p_val = stats.fisher_exact([[mutation_count - motif_mutation_count, motif_mutation_count], [ref_count - motif_count, motif_count]])
    p_value = p_val[1]
    chi_array = np.array([[mutation_count - motif_mutation_count, motif_mutation_count], [ref_count - motif_count, motif_count]])
    chi = stats.chi2_contingency(chi_array)[1]
    # if p_value <= 0.05:
    #  qvalues = multipletests(pvals=p_value, method='fdr_bh')
    #  print(qvalues)
    #  if qvalues[3] <= 0.05:
    #     print("significant")
    # print("odds_ratio: ", "p-value")
    return (p_val[0], p_value, chi)


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

    return zip(
        cycle([mutation[0]]),
        range(mutation[1] - range_size - motif_position, mutation[1] + len(motif) + range_size - motif_position),
        seq,
        cycle("+"))


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

    sample_list = defaultdict(list)

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
            sample_list[sample].append((chrom, pos, x, y))
            mutation = {
                'chrom': chrom,
                'pos': pos,
                'x': x,
                'y': y}
        return sample_list


def get_rev_comp_seq(sequence):
    # rev_comp_seq = "".join([complementary_nucleotide[i] for i in reversed(sequence)])
    rev_comp_seq = [(i[0], i[1], complementary_nucleotide[i[2]], "-") for i in reversed(sequence)]
    return rev_comp_seq


def mutated_base_forward(mutation, ref, alt):
    """
    :param mutation: [(record.CHROM, record.POS, record.REF, record.ALT)]
    :param ref: the nucleotide base pre-mutation
    :param alt: the nucleotide base post-mutation
    :return: how many mutations match specified ref and alt
    """
    assert ref != alt, "mutation should have different ref and alt nucleotides"
    assert len(ref) == 1 and len(alt) == 1, "ref and alt should be 1 nucleotide"

    comp_ref = complementary_nucleotide[ref]
    comp_alt = complementary_nucleotide[alt]
    # makes sure single base substitution
    if mutation[3] and mutation[2] and len(mutation[2]) == 1 \
            and len(mutation[3]) == 1 and len(mutation[3]) == 1 \
            and len(mutation[2]) == 1 and mutation[2] != mutation[3]:
        # mutation does not match the substitution
        if ref == mutation[2] and alt == mutation[3]:
            return True


def mutated_base_rev(mutation, ref, alt):
    """
    :param mutation: [(record.CHROM, record.POS, record.REF, record.ALT)]
    :param ref: the nucleotide base pre-mutation
    :param alt: the nucleotide base post-mutation
    :return: how many mutations match specified ref and alt
    """
    assert ref != alt, "mutation should have different ref and alt nucleotides"
    assert len(ref) == 1 and len(alt) == 1, "ref and alt should be 1 nucleotide"

    comp_ref = complementary_nucleotide[ref]
    comp_alt = complementary_nucleotide[alt]
    # makes sure single base substitution
    if mutation[3] and mutation[2] and len(mutation[2]) == 1 \
            and len(mutation[3]) == 1 and len(mutation[3]) == 1 \
            and len(mutation[2]) == 1 and mutation[2] != mutation[3]:
        # mutation does not match the substitution
        if comp_ref == mutation[2] and comp_alt == mutation[3]:
            return True


def find_matching_motifs(seq, motif):
    # print("Looking for motif {} in {}, {}".format(motif, sequence, len(sequence) - len(motif)))
    for i in range(len(seq) - len(motif) + 1):
        s = seq[i: i + len(motif)]
        for j, char in enumerate(motif):
            if s[j][2] not in bases_dict[char]:
                break
        else:
            yield seq[i]


def get_enrichment(filename, motif, motif_position, ref, alt, range_size, assembly=None):
    sample_list = read_MAF(filename)

    result_by_sample = defaultdict(list)

    assert range_size >= 0
    assert len(ref) == 1
    assert len(alt) == 1
    assert 0 <= motif_position < len(motif)

    # extra loop for sample in sample list
    for sample, mutation in sample_list.items():

        matching_bases = set()
        matching_motifs = set()
        matching_mutated_motifs = set()
        matching_mutated_bases = set()

        # extract the longest sequence we would ever need (motif + range_size)
        for mut in mutation:

            m = (mut[0], mut[1], mut[2])

            seq = list(motif_contexts(m, motif, motif_position, range_size, assembly))
            rev_seq = get_rev_comp_seq(seq)

            # not mutated:
            for data in seq:
                if data[2] == ref:
                    matching_bases.add(data)

            for motif_match in find_matching_motifs(seq, motif):
                matching_motifs.add(motif_match)

            # rev compl: not mutated:
            for item in rev_seq:
                if item[2] == ref:
                    matching_bases.add(item)

            for motif_match in find_matching_motifs(rev_seq, motif):
                matching_motifs.add(motif_match)

            # mutated:
            if mutated_base_forward(mut, ref, alt):
                m = (mut[0], mut[1], mut[2], "+")
                matching_mutated_bases.add(m)

                context_of_mutation = seq[range_size: range_size + len(motif)]
                for motif_match in find_matching_motifs(context_of_mutation, motif):
                    matching_mutated_motifs.add(motif_match)

            if mutated_base_rev(mut, ref, alt):
                m = (mut[0], mut[1], mut[2], "-")
                matching_mutated_bases.add(m)

                # rev compl:
                context_of_mutation = rev_seq[range_size: range_size + len(motif)]
                for motif_match in find_matching_motifs(context_of_mutation, motif):
                    matching_mutated_motifs.add(motif_match)

        motif_mutation_count = len(matching_mutated_motifs)
        mutation_count = len(matching_mutated_bases)
        ref_count = len(matching_bases)
        motif_count = len(matching_motifs)

        enrichment = (motif_mutation_count / mutation_count) / (motif_count / ref_count)

        if enrichment > 1:
            mut_load = (motif_mutation_count * (enrichment - 1)) / enrichment
        else:
            mut_load = 0.0
        p_val = get_stats(motif_mutation_count, mutation_count, motif_count, ref_count)
        result_by_sample[sample].append((round(mut_load), len(mutation)))

    return enrichment, round(mut_load), p_val[1], p_val[2], motif_mutation_count, \
           mutation_count, motif_count, ref_count, result_by_sample


def get_values(filename, motif, motif_position, ref, alt, assembly=None):
    mut_burden_count = []
    range_sizes = []
    with open("APOBEC_window.csv", "w") as csvfile: #create csv file
        filewriter = csv.writer(csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow([motif])
        filewriter.writerow(["window_size", "enrichment", "motif_mut_count", "mut_count", "motif_count", "ref_count", "mutation load", "fisher p-value", "chi p-value"])
        for val in range(0, 50 + 1, 2):
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
            filewriter.writerow([val*2 + len(motif), result, motif_mut, mut, motif_count, ref_count, mut_burden, pvalue, chi_pvalue])
        return list(zip(range_sizes, mut_burden_count))


def analyse_results_counts(filename, motif, motif_position, ref, alt, range_size, assembly=None):
    filer = get_enrichment(filename, motif, motif_position, ref, alt, range_size, assembly=assembly)[8]

    samples = []
    sig_mut = []
    non_sig_mut = []
    for key, val in filer.items():
        sig_mut.append(val[0][0])
        non_sig_mut.append(val[0][1] - val[0][0])
        samples.append(key)

    legend1 = "Number of Non-" + str(motif) + " Mutations"
    legend2 = "Number of " + str(motif) + " Mutations"

    p1 = plt.bar(range(len(filer)), non_sig_mut, align='center')
    p2 = plt.bar(range(len(filer)), sig_mut, bottom = non_sig_mut, align='center', color = "crimson")

    plt.margins(0.2)
    plt.xlabel("Sample ID")
    plt.ylabel('Mutational Burden')
    plt.title('Mutational Burden By Sample')
    plt.xticks(range(len(filer)), filer.keys(), rotation = '45', ha = "right")
    plt.legend((p1[0], p2[0]), (legend1, legend2))
    plt.savefig("mutation_attribution.jpg", bbox_inches = "tight")
    plt.show()
    return filer


def analyse_results_prop(filename, motif, motif_position, ref, alt, range_size, assembly=None):
    filer = get_enrichment(filename, motif, motif_position, ref, alt, range_size, assembly=assembly)[8]

    samples = []
    percent_sig = []
    non_sig = []
    for key, val in filer.items():
        percent_sig.append(val[0][0]/val[0][1])
        non_sig.append(1-val[0][0]/val[0][1])
        samples.append(key)

    legend2 = "Percentage of " + str(motif) + " Mutations"
    legend1 = "Percentage of Non-" + str(motif) + " Mutations"

    p1 = plt.bar(range(len(filer)), non_sig, align='center')
    p2 = plt.bar(range(len(filer)), percent_sig, bottom = non_sig, align='center', color="crimson")

    plt.margins(0.2)
    plt.xlabel("Sample ID")
    plt.ylabel('Percentage of Matching Mutations')
    plt.title('Percentage of Matching Mutations by Sample')
    plt.xticks(range(len(filer)), filer.keys(), rotation='45', ha = "right")
    plt.legend((p1[0], p2[0]), (legend1, legend2))
    #axes = plt.gca()
    plt.savefig("mutation_prop.jpg", bbox_inches="tight")
    plt.show()
    return filer


def graph_range():
    x = [2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 74, 78 ,82, 86, 90, 94, 98, 102]
    y = [3,
            2,
            5,
            5,
            6,
            7,
            7,
            7,
            7,
            7,
            7,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8]
    plt.plot(x, y, 'ro')
    for i in np.arange(0, len(x), 1):
        plt.plot(x[i:i + 2], y[i:i + 2], 'k-')
    plt.axis([0, 105, 0, 10])
    plt.title("Mutational Burden Depends Upon Window Size")
    plt.xlabel("Window Size")
    plt.ylabel("Mutational Burden")
    plt.savefig("rangesize1_graph.png")
    plt.show(block=True)
    return "graph completed"

def graph_motif_sig():
    samples = ["AGEING Signature", "AGEING Motif"]
    values = [22, 11]
    p1 = plt.bar(samples, values, align='center', color = 'crimson')
    plt.margins(0.2)
    plt.xlabel("Analysis Method")
    plt.ylabel('Mutational Burden')
    plt.title('Mutational Burden Depends on Analysis Method')
    #axes = plt.gca()
    plt.savefig("motif_sig_cpg.jpg", bbox_inches="tight")
    plt.show()

def get_rev_comp(sequence):
    """
    :param sequence: DNA forward sequence
    return: reverse complementary DNA sequence
     """
    rev_comp_seq = ""
    for i in reversed(sequence):
        rev_comp_seq += complementary_nucleotide[i]
    return rev_comp_seq

if __name__ == '__main__':
    #pprint(graph_motif_sig())
    #pprint(get_values("tcga_A0C8.maf", "TCW", 1, "C", "T", assembly=37))
    #pprint(get_values("tcga_A0C8.maf", "TCW", 1, "C", "G", assembly=37))
    print(get_rev_comp("AATGCTAGCTAGCTAGCTAGCTAGCTGATGCTAGCTAGCTAGCTGATCGT"))