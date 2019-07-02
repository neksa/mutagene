import vcf
import twobitreader as tbr
from collections import defaultdict
from statsmodels.stats.multitest import multipletests

import scipy.stats as stats
import csv
import numpy as np
from itertools import cycle


nucleotides = "ACGT"  # this order of nucleotides is important for reversing
complementary_nucleotide = dict(zip(nucleotides, reversed(nucleotides)))
complementary_nucleotide['N'] = 'N'

TWOBIT_GENOMES_PATH = '/net/pan1/mutagene/data/genomes/'

bases_dict = {"A": "A", "G": "G", "T": "T", "C": "C", "W": "AT", "S": "CG", "M": "AC", "K": "GT", "R": "AG", "Y": "CT",
              "B": "TCG", "D": "AGT", "H": "ACT", "V": "ACG", "N": "AGTC"}

comp_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "W": "AT", "S": "CG", "K": "AC", "M": "GT", "Y": "AG",
             "R": "CT", "V": "TCG", "H": "AGT", "D": "ACT", "B": "ACG", "N": "ATGC"}


def get_stats(motif_mutation_count, mutation_count, motif_count, ref_count):
    """
    :param motif_mutation_count: number of mutations that match the motif
    :param mutation_count: number of mutations that match the specified ref and alt nucleotides
    :param motif_count: number of motif occurrences in seq
    :param ref_count: number of reference nucleotide occurrences in seq
    :return: p-value and odds ratio
    """
    p_val = stats.fisher_exact([[mutation_count - motif_mutation_count, motif_mutation_count], [ref_count - motif_count, motif_count]])
    p_value = p_val[1]
    chi_array = np.array([[mutation_count - motif_mutation_count, motif_mutation_count], [ref_count - motif_count, motif_count]])
    chi = stats.chi2_contingency(chi_array)[1]
    # if p_value <= 0.05:
    #  qvalues = multipletests(pvals=p_value, method='fdr_bh')
    #  print(qvalues)
    #  if qvalues[3] <= 0.05:
    #     print("significant")
    #print("odds_ratio: ", "fisher p-value", "chi p-value" )
    return (p_val[0], p_value, chi)


def motif_contexts(mutation, motif, motif_position, range_size, assembly=None):
    """

    :param mutations:
    :param motif: contains the ref nucleotide and surrounding nucleotides
    :param range_size: contains the ref nucleotide and surrounding nucleotides
    :param assembly: specific genome sequence, default 38
    :return: number of times a motif has occured in a DNA sequence, incl. reverse complementary DNA strand
    """
    contexts = {}
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

    if assembly not in twobit_files:
        return contexts
    else:
        twobit_file = genomes_path + "/" + twobit_files[assembly] + ".2bit"
        f = tbr.TwoBitFile(twobit_file)

        start = mutation[1] - 1  # zero-based numbering; start=0
        chrom = str(mutation[0])
        chromosome = chrom if chrom.startswith('chr') else 'chr' + chrom
        chromosome = chromosome_name_mapping.get(chromosome, chromosome)

        if chromosome in f:
            try:
                seq = f[chromosome][start - range_size - motif_position:
                                    start + len(motif) + range_size - motif_position]

                seq = seq.upper()

            except Exception as e:
                print("TwoBit exception", str(e))
        else:
            print("NO CHROM", chromosome)
            pass

    return list(zip(
        cycle([mutation[0]]),
        range(mutation[1] - range_size - motif_position, mutation[1] + len(motif) + range_size - motif_position),
        seq,
        cycle("+")))


def read_VCF(filename):
    """
    :param filename: name of VCF file in quotes
    :return: lists of single base substitution mutations in VCF
    """
    vcf_reader = vcf.Reader(filename=filename)

    mutations_samples = defaultdict(list)

    for record in vcf_reader:
        mutations_samples[record.ID].append((record.CHROM, record.POS, record.REF, record.ALT))
    return mutations_samples


def get_rev_comp_seq(sequence):
    """
    :param sequence: DNA forward sequence
    return: reverse complementary DNA sequence
     """
    rev_comp_seq = [(i[0], i[1], complementary_nucleotide[i[2]], "-") for i in reversed(sequence)]

    return rev_comp_seq


def mutated_base_forward(mutations, ref, alt):

    assert ref != alt, "mutation should have different ref and alt nucleotides"

    assert len(ref) == 1 and len(alt) == 1, "ref and alt should be 1 nucleotide"

    # make sure it is a single base substitution
    if mutations[3][0] and mutations[2] and len(mutations[2]) == 1 \
            and len(mutations[3][0]) == 1 and len(mutations[3]) == 1 \
            and len(mutations[2]) == 1 and mutations[2] != mutations[3][0]:

        if ref == mutations[2] and alt == mutations[3][0]:
            return True


def mutated_base_rev(mutations, ref, alt):

    assert ref != alt, "mutation should have different ref and alt nucleotides"

    assert len(ref) == 1 and len(alt) == 1, "ref and alt should be 1 nucleotide"

    comp_ref = complementary_nucleotide[ref]
    comp_alt = complementary_nucleotide[alt]
    # make sure it is a single base substitution
    if mutations[3][0] and mutations[2] and len(mutations[2]) == 1 \
            and len(mutations[3][0]) == 1 and len(mutations[3]) == 1 \
            and len(mutations[2]) == 1 and mutations[2] != mutations[3][0]:

        if comp_ref == mutations[2] and comp_alt == mutations[3][0]:
            return True


def find_matching_motifs(seq, motif, motif_position):
    # print("Looking for motif {} in {}, {}".format(motif, sequence, len(sequence) - len(motif)))
    for i in range(len(seq) - len(motif) + 1):
        s = seq[i: i + len(motif)]
        for j, char in enumerate(motif):
            if s[j][2] not in bases_dict[char]:
                break
        else:
            yield seq[i + motif_position]


def get_enrichment(filename, motif, motif_position, ref, alt, range_size, assembly=None):
    """
    :param filename: name of VCF file in quotes
    :param motif: mutated base plus nucleotide context
    :param motif_position: zero-based numbering, which nucleotide is being mutated in motif
    :param ref: the nucleotide base pre-mutation
    :param alt: the nucleotide base post-mutation
    :param range_size: the number of nucleotides the sequence extends to excluding the motif on each side
    :param assembly: specific genome sequence, default 38
    :return: enrichment and mutation load
        """
    mutation_samples = read_VCF(filename)
    result_by_sample = defaultdict(list)

    assert range_size >= 0
    assert len(ref) == 1
    assert len(alt) == 1
    assert 0 <= motif_position < len(motif)


    matching_bases = set()
    matching_motifs = set()
    matching_mutated_motifs = set()
    matching_mutated_bases = set()

    #extra loop for sample in sample list
    for sample, mutation in mutation_samples.items():

        for mut in mutation:
            # extract the longest sequence we would ever need (motif + range_size)
            m = (mut[0], mut[1], mut[2])
            seq = motif_contexts(m, motif, motif_position, range_size, assembly)

            rev_seq = get_rev_comp_seq(seq)

            # print("--------------------")
            # pprint(mutation)
            # print(seq)
            # print(rev_seq)

            #not mutated
            for data in seq:
                if data[2] in bases_dict[ref]:
                    matching_bases.add(data)

            for motif_match in find_matching_motifs(seq, motif, motif_position):
                matching_motifs.add(motif_match)

            # rev compl: not mutated:
            for item in rev_seq:
                if item[2] in bases_dict[ref]:
                    matching_bases.add(item)

            for motif_match in find_matching_motifs(rev_seq, motif, len(motif) -1 - motif_position):
                matching_motifs.add(motif_match)

            # mutated:
            if mutated_base_forward(mut, ref, alt):
                m = (mut[0], mut[1], mut[2], "+")
                matching_mutated_bases.add(m)

                context_of_mutation = seq[range_size: range_size + len(motif)]
                for motif_match in find_matching_motifs(context_of_mutation, motif, motif_position):
                    matching_mutated_motifs.add(motif_match)

            if mutated_base_rev(mut, ref, alt):
                m = (mut[0], mut[1], mut[2], "-")
                matching_mutated_bases.add(m)

                # rev compl:
                context_of_mutation = rev_seq[range_size: range_size + len(motif)]
                for motif_match in find_matching_motifs(context_of_mutation, motif, len(motif) -1 - motif_position):
                    matching_mutated_motifs.add(motif_match)

        motif_mutation_count = len(matching_mutated_motifs)
        mutation_count = len(matching_mutated_bases)
        ref_count = len(matching_bases)
        motif_count = len(matching_motifs)

        try:
            enrichment = (motif_mutation_count / mutation_count) / (motif_count / ref_count)
            if enrichment > 1:
                mut_load = (motif_mutation_count * (enrichment - 1)) / (enrichment)
            else:
                mut_load = 0.0
            p_val = get_stats(motif_mutation_count, mutation_count, motif_count, ref_count)
            result_by_sample[sample].append((round(mut_load)))
        except:
            "values too low to calc enrichment"
            mut_load = 0.0
    try:
        return enrichment, round(mut_load), p_val[1], p_val[2], motif_mutation_count, mutation_count, motif_count, ref_count, (round(mut_load))/len(mutation_samples), result_by_sample
    except:
        return 0.0


def get_values(filename, motif, motif_position, ref, alt, assembly=None):
    range_sizes = []
    with open("APOBEC_window1.csv", "w") as csvfile: #create csv file
        filewriter = csv.writer(csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow([motif])
        filewriter.writerow(["range_size", "enrichment", "motif_mut_count", "mut_count", "motif_count", "ref_count", "mutation load", "fisher p-value", "chi p-value"])
        for val in range(0, 5000 + 1, 100):
            range_sizes.append(val * 2 + len(motif))
            analysis = get_enrichment(filename, motif, motif_position, ref, alt, val, assembly=assembly)
            result = analysis[0]
            mut_burden = analysis[1]
            pvalue = analysis[2]
            chi_pvalue = analysis[3]
            motif_mut = analysis[4]
            mut = analysis[5]
            motif_count = analysis[6]
            ref_count = analysis[7]
            filewriter.writerow([val * 2 + len(motif), result, motif_mut, mut, motif_count, ref_count, mut_burden, pvalue, chi_pvalue])
        return "completed"


if __name__ == '__main__':
    #print(get_enrichment("data.vcf", "TCW", 1, "C", "T", 29, assembly=37))
    print(get_enrichment("data.vcf", "CG", 0, "C", "T", 19, assembly=38))

    #print(get_enrichment("data.vcf", "TCG", 1, "C", "T", 29, assembly=37))
    #print(special_motifs_enrichment('data.vcf','uv', assembly=37))

# add adaptor to put mutated base in caps
