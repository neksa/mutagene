
import json
from subprocess import Popen, PIPE

try:
    from .identify_signature import *
    # from . import *
    from .mutations import *
    from .json_structs import *
except:
    pass

from itertools import chain, islice

from collections import namedtuple, defaultdict
import gzip
import csv
import os
# from .COSMIC import iter_CLP_mutations


TWOBIT_GENOMES_PATH = '/net/pan1/mutagene/data/genomes/'


def make_NCI60_report(prefix, n=30, decomposition='L'):
    s = {5: 'MUTAGENE A.', 10: 'MUTAGENE B.', 30: 'COSMIC '}[n]

    COSMIC_CLINICAL = prefix + "/cellline_nci60_clinical_data.tsv"

    stats_columns = ("N", "Processed_C", "Processed_NC", "N_C", "N_NC", "Skipped_C", "Skipped_NC", "Skipped_indels_C", "Skipped_indels_NC", "Skipped_nucleotide_C", "Skipped_nucleotide_NC", "Skipped_context_C", "Skipped_context_NC")
    with open(COSMIC_CLINICAL) as f, open(prefix + "/report_{}{}.txt".format(decomposition, n), 'w') as o:
        o.write("SAMPLE_ID\tCANCER_TYPE\t" + "\t".join(stats_columns) + "\t" + "\t".join(["{}{}".format(s, x) for x in range(1, n + 1)]) + "\tOther\tUnexplained\n")
        for line in islice(f, 1, None):
            if line.startswith('#'):
                continue
            if line.startswith('SAMPLE_ID'):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            SAMPLE_ID = fields[0]
            CANCER = fields[3]

            if not os.path.isfile(prefix + "/profiles/{}.stats".format(SAMPLE_ID)):
                continue

            print(SAMPLE_ID, CANCER)

            o.write("{}\t{}\t".format(SAMPLE_ID, CANCER))

            stats_values = {}
            with open(prefix + "/profiles/{}.stats".format(SAMPLE_ID)) as stats:
                for line in stats:
                    stats, value = line.split("\t")
                    stats_values[stats] = value.strip()

            for stats in stats_columns:
                o.write("{}\t".format(stats_values[stats]))

            exposure = defaultdict(float)
            with open(prefix + "/decompose/{}_{}-{}.txt".format(SAMPLE_ID, decomposition, n)) as d:
                for line in d:
                    signature, value = line.split("\t")
                    value = float(value)
                    if signature.startswith(s):
                        exposure[signature] = value
                    if signature.startswith("Other"):
                        exposure['other'] = value
            unexplained = 1.0 - sum(exposure.values())
            o.write("\t".join(["{0:.2f}".format(exposure["{}{}".format(s, x)]) for x in range(1, n + 1)]) + "\t{:.2f}\t{:.2f}\n".format(exposure['other'], unexplained))


def read_MAF_extended(muts, asm=None):
    cn = complementary_nucleotide
    N_skipped = 0

    samples = defaultdict(list)

    raw_mutations = []
    for i, line in enumerate(muts.split("\n")):
        if line.startswith("#"):
            continue
        if line.startswith("Hugo_Symbol"):
            continue
        if len(line) < 10:
            continue

        col_list = line.split("\t")
        # for a, b in enumerate(col_list):
        #     print(a, b)

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

            pos = int(col_list[5])      # MAF POS START
            pos_end = int(col_list[6])  # MAF POS END
            x = col_list[10]            # MAF REF
            y1 = col_list[11]           # MAF ALT1
            y2 = col_list[12]           # MAF ALT2
        except:
            # raise
            continue

        sample = None
        try:
            sample = col_list[15]           # Tumor_Sample_Barcode
        except:
            raise

        if pos != pos_end:
            continue

        # skip if found unexpected nucleotide characters
        if len(set([x, y1, y2]) - set(nucleotides)) > 0:
            continue

        y = y1 if y1 != x else None
        y = y2 if y2 != x else y
        if y is None:
            continue

        samples[sample].append((chrom, pos, x, y))

    processing_stats = defaultdict(dict)

    MAX = 100000000
    for sample, raw_mutations in samples.items():
        print("sample", sample, len(raw_mutations))

        # if sample != "786_0":
        #     break

        mutations = defaultdict(float)
        if len(raw_mutations) > 0:
            if len(raw_mutations) > MAX:
                raw_mutations = raw_mutations[:MAX]

            contexts = get_context_batch(raw_mutations, asm)

            if contexts is None:
                return None, None

            if len(contexts) == 0:
                return None, None

            for (chrom, pos, x, y) in raw_mutations:
                p5, p3 = contexts.get((chrom, pos), ("N", "N"))

                if len(set([p5, x, y, p3]) - set(nucleotides)) > 0:
                    # print("Skipping invalid nucleotides")
                    N_skipped += 1
                    continue

                if x in "CT":
                    mutations[p5 + p3 + x + y] += 1.0
                else:
                    # complementary mutation
                    mutations[cn[p3] + cn[p5] + cn[x] + cn[y]] += 1.0

        N_loaded = int(sum(mutations.values()))
        processing_stats[sample] = {'loaded': N_loaded, 'skipped': N_skipped, 'format': 'MAF'}
        samples[sample] = mutations

    return samples, processing_stats


uniq_mut = defaultdict(set)


def read_CLP_extended(f, asm=None, ids=None):
    cn = complementary_nucleotide
    # N_skipped = 0
    # uniq_mut = defaultdict(set)

    samples = defaultdict(list)
    samples_raw = defaultdict(list)
    samples_mutations = defaultdict(int)
    samples_skipped_indels = defaultdict(int)
    samples_skipped_nucleotide = defaultdict(int)
    samples_skipped_context = defaultdict(int)

    CLPRecord = namedtuple('CLPRecord',
        'gene, accession, CDS_length, HGNC, sample_name, sample_id, tumor_id, site, site_subtype1, site_subtype2, site_subtype3, ' +
        'histology, histology_subtype1, histology_subtype2, histology_subtype3, GWS, mut_id, mut_CDS, mut_AA, mut_descr, mut_zygosity, LOH, ' +
        'GRCh, GRC_pos, GRC_strand, SNP, FATHMM, FATHMM_score, somatic, verification, PMID, study_id, institute, institute_address, catalogue_no, ' +
        'sample_source, tumor_origin, age')

    raw_mutations = []
    for rec in map(CLPRecord._make, csv.reader(islice(f, 1, None), delimiter='\t')):
        # print rec.gene
        # print rec.mut_CDS
        # print rec.mut_AA
        # print rec.somatic

        try:
            sample = rec.sample_name.replace("-", "_").upper()

            chrom = ''
            pos = 0
            pos_end = 0
            if len(rec.GRC_pos) > 0:
                chrom, start_end = rec.GRC_pos.split(":")
                pos, pos_end = start_end.split("-")
                pos = int(pos)
                pos_end = int(pos_end)

            x = rec.mut_CDS[-3]
            y = rec.mut_CDS[-1]

            if ids and sample not in ids:
                continue

            mut = (rec.GRC_pos if len(rec.GRC_pos) > 0 else rec.mut_id, x, y)
            if mut in uniq_mut[sample]:
                # It's probably another isoform!
                # print("SKIPPING DUPLICATE", sample, rec.GRC_pos, x, y)
                continue
            uniq_mut[sample] |= set([mut])

            samples_mutations[sample] += 1

        except Exception as e:
            # raise
            print("SKIPPING RECORD", str(e), rec)
            continue

        if pos != pos_end or "ins" in rec.mut_CDS or "del" in rec.mut_CDS:
            samples_skipped_indels[sample] += 1
            continue

        if len(chrom) == 0:
            print("SKIPPING. NO GRCh specified", rec)
            samples_skipped_nucleotide[sample] += 1
            continue

        # skip if found unexpected nucleotide characters
        if len(set([x, y]) - set(nucleotides)) > 0:
            samples_skipped_nucleotide[sample] += 1
            continue

        samples[sample].append((chrom, pos, x, y))

    processing_stats = defaultdict(dict)

    # MAX = 100000000
    for sample, raw_mutations in samples.items():
        print("sample", sample, len(raw_mutations))

        # if sample != "786_0":
        #     break

        mutations = defaultdict(float)
        if len(raw_mutations) > 0:
            # if len(raw_mutations) > MAX:
            #     raw_mutations = raw_mutations[:MAX]

            contexts = get_context_batch(raw_mutations, asm)

            if contexts is None:
                return None, None

            if len(contexts) == 0:
                return None, None

            for (chrom, pos, x, y) in raw_mutations:
                p5, p3 = contexts.get((chrom, pos), ("N", "N"))

                if len(set([p5, x, y, p3]) - set(nucleotides)) > 0:
                    print("Skipping invalid nucleotides {}:{} {}[{}>{}]{}".format(chrom, pos, p5, x, y, p3))
                    # N_skipped += 1
                    samples_skipped_context[sample] += 1
                    continue

                if x in "CT":
                    mutations[p5 + p3 + x + y] += 1.0
                else:
                    # complementary mutation
                    mutations[cn[p3] + cn[p5] + cn[x] + cn[y]] += 1.0

        N = samples_mutations[sample]
        N_loaded = int(sum(mutations.values()))
        processing_stats[sample] = {'processed': N, 'loaded': N_loaded, 'skipped': N - N_loaded, 'skipped_indels': samples_skipped_indels[sample], 'skipped_context': samples_skipped_context[sample], 'skipped_nucleotide': samples_skipped_nucleotide[sample] }
        samples[sample] = mutations
        samples_raw[sample] = raw_mutations

    return samples, processing_stats, samples_raw


def read_CLP_NCV(f, asm=None, ids=None):
    cn = complementary_nucleotide
    # N_skipped = 0

    # uniq_mut = defaultdict(set)
    samples = defaultdict(list)
    samples_raw = defaultdict(list)
    samples_mutations = defaultdict(int)
    samples_skipped_nucleotide = defaultdict(int)
    samples_skipped_context = defaultdict(int)
    samples_skipped_indels = defaultdict(int)

    CLPNCVRecord = namedtuple('CLPNCVRecord',
        'sample_name, sample_id, tumor_id, site, site_subtype1, site_subtype2, site_subtype3, ' +
        'histology, histology_subtype1, histology_subtype2, histology_subtype3, GWS, mut_id, mut_zygosity, ' +
        'GRCh, GRC_pos, somatic, WT, MUT, SNP, FATHMM_score, FATHMM_groups, FATHMM_MKL_score, FATHMM_MKL_groups, genome, exome, study')

    raw_mutations = []
    for rec in map(CLPNCVRecord._make, csv.reader(islice(f, 1, None), delimiter='\t')):
        try:
            sample = rec.sample_name.replace("-", "_").upper()

            chrom, start_end = rec.GRC_pos.split(":")
            pos, pos_end = start_end.split("-")
            pos = int(pos)
            pos_end = int(pos_end)

            x = rec.WT
            y = rec.MUT

            if ids and sample not in ids:
                continue

            samples_mutations[sample] += 1

            mut = (rec.GRC_pos, x, y)
            if mut in uniq_mut[sample]:
                # It's probably another isoform!
                print("SKIPPING DUPLICATE", sample, rec.GRC_pos, x, y)
                continue
            uniq_mut[sample] |= set([mut])

        except Exception as e:
            # raise
            print("SKIPPING RECORD", str(e), rec)
            continue

        if pos != pos_end or len(x) != 1 or len(y) != 1:
            # print("SKIPPING INDEL", pos, pos_end, x, y)
            samples_skipped_indels[sample] += 1
            continue

        # skip if found unexpected nucleotide characters
        if len(set([x, y]) - set(nucleotides)) > 0:
            # print("SKIPPING UNEXPECTED", x, y)
            samples_skipped_nucleotide[sample] += 1
            continue

        samples[sample].append((chrom, pos, x, y))

    processing_stats = defaultdict(dict)

    # MAX = 100000000
    for sample, raw_mutations in samples.items():
        print("sample", sample, len(raw_mutations))

        # if sample != "786_0":
        #     break

        mutations = defaultdict(float)
        if len(raw_mutations) > 0:
            # if len(raw_mutations) > MAX:
            #     raw_mutations = raw_mutations[:MAX]

            contexts = get_context_batch(raw_mutations, asm)

            if contexts is None:
                return None, None

            if len(contexts) == 0:
                return None, None

            for (chrom, pos, x, y) in raw_mutations:
                p5, p3 = contexts.get((chrom, pos), ("N", "N"))

                if len(set([p5, x, y, p3]) - set(nucleotides)) > 0:
                    # print("Skipping invalid nucleotides")
                    print("Skipping invalid nucleotides {}[{}>{}]{}".format(p5, x, y, p3))
                    # N_skipped += 1
                    samples_skipped_context[sample] += 1
                    continue

                if x in "CT":
                    mutations[p5 + p3 + x + y] += 1.0
                else:
                    # complementary mutation
                    mutations[cn[p3] + cn[p5] + cn[x] + cn[y]] += 1.0

        N = samples_mutations[sample]
        N_loaded = int(sum(mutations.values()))
        # processing_stats[sample] = {'loaded': N_loaded, 'skipped': N - N_loaded}
        processing_stats[sample] = {'processed': N, 'loaded': N_loaded, 'skipped': N - N_loaded, 'skipped_indels': samples_skipped_indels[sample], 'skipped_context': samples_skipped_context[sample], 'skipped_nucleotide': samples_skipped_nucleotide[sample] }
        samples[sample] = mutations
        samples_raw[sample] = mutations_raw

    return samples, processing_stats, samples_raw



def deconstruct_sigs(profile_fname, sample):
    script = """
library(jsonlite)
library(deconstructSigs)
s <- t(read.table('{}', sep="\t", header=FALSE, row.names=1))
row.names(s) <- '{}'
s <- as.data.frame(s)
w <- whichSignatures(s/sum(s), signatures.ref=signatures.cosmic)
toJSON(w)
"""
    script = script.format(profile_fname, sample).encode("utf-8")
    proc = Popen(["/Users/gonceare/anaconda3/envs/mutagene/bin/Rscript", "-"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    out, err = proc.communicate(script)
    # exitcode = proc.returncode
    # print(exitcode, out, err)
    json_string = out.decode("utf-8")
    w = json.loads(json_string)
    # pprint.pprint(w)
    result = []
    for k, v in w['weights'][0].items():
        if k.startswith('_row'):
            continue
        if float(v) == 0.0:
            continue
        result.append({
            'name': k.replace('Signature.', 'COSMIC '),
            'score': v})
    return None, result

# GRC37
# from cosmic website 468 noncoding + 330 coding = 798
# grep files: 519 + 350
# 818 loaded 51 skipped = 869

# GRC38
# ??? 394
# 321 + 461 = 782
# grep 344 + 511 = 855
#OLD ISO: BT_549': {'loaded': 805, 'skipped': 50} = 855
# 'BT_549': {'loaded': 714, 'skipped': 36}: = 750
#


def export_maf(prefix, maf_file):
    """Create a MAF file with mutations

    Mutations in NCI-60 cell lines are extracted from COSMIC CLP

    Arguments:
        prefix {string} -- Location of COSMIC data
        maf_file {string} -- output file name

    Returns:
        None
    """

    """
    GDC MAF Format v.1.0.0
    From https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/

    1 - Hugo_Symbol   HUGO symbol for the gene (HUGO symbols are always in all caps). "Unknown" is used for regions that do not correspond to a gene
    2 - Entrez_Gene_Id    Entrez gene ID (an integer). "0" is used for regions that do not correspond to a gene region or Ensembl ID
    3 - Center  One or more genome sequencing center reporting the variant
    4 - NCBI_Build  The reference genome used for the alignment (GRCh38)
    5 - Chromosome  The affected chromosome (chr1)
    6 - Start_Position  Lowest numeric position of the reported variant on the genomic reference sequence. Mutation start coordinate
    7 - End_Position    Highest numeric genomic position of the reported variant on the genomic reference sequence. Mutation end coordinate
    8 - Strand  Genomic strand of the reported allele. Currently, all variants will report the positive strand: '+'
    9 - Variant_Classification  Translational effect of variant allele
    10 - Variant_Type   Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP (di-nucleotide polymorphism) but for three consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of four or more (SNP, DNP, TNP, ONP, INS, DEL, or Consolidated)
    11 - Reference_Allele   The plus strand reference allele at this position. Includes the deleted sequence for a deletion or "-" for an insertion
    12 - Tumor_Seq_Allele1  Primary data genotype for tumor sequencing (discovery) allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases
    13 - Tumor_Seq_Allele2  Tumor sequencing (discovery) allele 2
    14 - dbSNP_RS   The rs-IDs from the   dbSNP database, "novel" if not found in any database used, or null if there is no dbSNP record, but it is found in other databases
    15 - dbSNP_Val_Status   The dbSNP validation status is reported as a semicolon-separated list of statuses. The union of all rs-IDs is taken when there are multiple
    16 - Tumor_Sample_Barcode   Aliquot barcode for the tumor sample
    17 - Matched_Norm_Sample_Barcode    Aliquot barcode for the matched normal sample
    """

    COSMIC_CLP = prefix + "/38/CosmicCLP_MutantExport.tsv.gz"
    COSMIC_CLP_NCV = prefix + "/38/CosmicCLP_NCVExport.tsv.gz"
    COSMIC_CLINICAL = prefix + "/cellline_nci60_clinical_data.tsv"

    ids = []
    with open(COSMIC_CLINICAL) as f:
        for line in islice(f, 1, None):
            fields = line.split("\t")
            ids.append(fields[0])

    with open(maf_file, 'w') as o:
        with gzip.open(COSMIC_CLP, mode='rt', encoding='UTF-8') as f:
            _, _, samples_raw = read_CLP_extended(f, 38, ids)
        with gzip.open(COSMIC_CLP_NCV, mode='rt', encoding='UTF-8') as f:
            _, _, samples_raw_NCV = read_CLP_NCV(f, 38, ids)

        for sample in samples_raw.keys():
            for mutation in chain(samples_raw[sample], samples_raw_NCV[sample]):
                chrom, pos, x, y = mutation
                o.write("\t\t\tGRCh38\t{}\t{}\t{}\t+\t\tSNP\t{}\t{}\t{}\t\t\t{}\t{}\n".format(
                    chrom, pos, pos, x, y, y, sample, sample))


def analyze_nci60_samples(prefix):
    if False:
        MAF_file = prefix + '/data_mutations_extended.txt'
        with open(MAF_file) as f:
            MAF = f.read()
            # print(MAF)
            samples, processing_stats = read_MAF_extended(MAF, 37)

    COSMIC_CLP = prefix + "/38/CosmicCLP_MutantExport.tsv.gz"
    COSMIC_CLP_NCV = prefix + "/38/CosmicCLP_NCVExport.tsv.gz"
    COSMIC_CLINICAL = prefix + "/cellline_nci60_clinical_data.tsv"
    THRESHOLD = 0.0  # 0.05

    ids = []
    with open(COSMIC_CLINICAL) as f:
        for line in islice(f, 1, None):
            fields = line.split("\t")
            ids.append(fields[0])

    with gzip.open(COSMIC_CLP, mode='rt', encoding='UTF-8') as f:
        samples, processing_stats = read_CLP_extended(f, 38, ids)
    # print(processing_stats)

    with gzip.open(COSMIC_CLP_NCV, mode='rt', encoding='UTF-8') as f:
        samples_NCV, processing_stats_NCV = read_CLP_NCV(f, 38, ids)
    # print(processing_stats_NCV)

    for k in samples_NCV.keys():
        for i in samples_NCV[k].keys():
            samples[k][i] += samples_NCV[k][i]

    # for k in processing_stats_NCV.keys():
    #     processing_stats[k]['loaded'] += processing_stats_NCV[k]['loaded']
    #     processing_stats[k]['skipped'] += processing_stats_NCV[k]['skipped']

    def format_numbers(d):
        d['score'] = "{:6.4g}".format((d['score']))
        return d

    # print(samples)
    # return

    for sample, mutations in samples.items():
        SAMPLE_ID = str(sample)
        # if sample != "786_0":
        #     break

        mutational_profile = get_mutational_profile(mutations)
        mutational_profile_counts = get_mutational_profile(mutations, counts=True)
        query_signature = make_bar_struct_from_values(mutational_profile)
        query_formatted = format_profile(mutational_profile)

        oname = prefix + "/profiles/{}.profile".format(SAMPLE_ID)
        print("ONAME:", oname)
        with open(oname, 'w') as o:
            query_formatted = format_profile(mutational_profile)
            o.write(query_formatted)

        oname = prefix + "/profiles/{}.counts".format(SAMPLE_ID)
        with open(oname, 'w') as o:
            query_formatted = format_profile(mutational_profile_counts, counts=True)
            o.write(query_formatted)

        oname = prefix + "/profiles/{}.stats".format(SAMPLE_ID)
        # print("ONAME:", oname)
        with open(oname, 'w') as o:
            o.write("N\t{}\n".format(int(sum(mutational_profile_counts))))
            o.write("Processed_C\t{}\n".format(processing_stats[sample]['processed']))
            o.write("Processed_NC\t{}\n".format(processing_stats_NCV[sample]['processed']))
            o.write("N_C\t{}\n".format(processing_stats[sample]['loaded']))
            o.write("N_NC\t{}\n".format(processing_stats_NCV[sample]['loaded']))
            o.write("Skipped_C\t{}\n".format(processing_stats[sample]['skipped']))
            o.write("Skipped_NC\t{}\n".format(processing_stats_NCV[sample]['skipped']))
            o.write("Skipped_indels_C\t{}\n".format(processing_stats[sample]['skipped_indels']))
            o.write("Skipped_indels_NC\t{}\n".format(processing_stats_NCV[sample]['skipped_indels']))
            o.write("Skipped_nucleotide_C\t{}\n".format(processing_stats[sample]['skipped_nucleotide']))
            o.write("Skipped_nucleotide_NC\t{}\n".format(processing_stats_NCV[sample]['skipped_nucleotide']))
            o.write("Skipped_context_C\t{}\n".format(processing_stats[sample]['skipped_context']))
            o.write("Skipped_context_NC\t{}\n".format(processing_stats_NCV[sample]['skipped_context']))

        classification_method = "rf"
        cancer_type_matches = [format_numbers(d) for d in identify_cancer_type(mutational_profile, classification_method)]
        primary_site_matches = [format_numbers(d) for d in identify_primary_site(mutational_profile, classification_method)]

        for method in "LJRZ":
        # for method in "D":
            # if method == 'M':
            #     _, contributing_signatures_A = decompose_mutational_profile(mutational_profile, "CLA")
            #     _, contributing_signatures_B = decompose_mutational_profile(mutational_profile, "CLB")
            #     _, contributing_signatures_cosmic = decompose_mutational_profile(mutational_profile, "IS")
            if method == 'L':
                _, _, contributing_signatures_A = decompose_mutational_profile_counts(mutational_profile_counts, "CLA", 'MLE', others_threshold=THRESHOLD)
                _, _, contributing_signatures_B = decompose_mutational_profile_counts(mutational_profile_counts, "CLB", 'MLE', others_threshold=THRESHOLD)
                _, _, contributing_signatures_cosmic = decompose_mutational_profile_counts(mutational_profile_counts, "IS", 'MLE', others_threshold=THRESHOLD)
            if method == 'J':
                _, _, contributing_signatures_A = decompose_mutational_profile_counts(mutational_profile_counts, "CLA", 'js')
                _, _, contributing_signatures_B = decompose_mutational_profile_counts(mutational_profile_counts, "CLB", 'js')
                _, _, contributing_signatures_cosmic = decompose_mutational_profile_counts(mutational_profile_counts, "IS", 'js')
            if method == 'R':
                _, _, contributing_signatures_A = decompose_mutational_profile_counts(mutational_profile_counts, "CLA", 'frobenius')
                _, _, contributing_signatures_B = decompose_mutational_profile_counts(mutational_profile_counts, "CLB", 'frobenius')
                _, _, contributing_signatures_cosmic = decompose_mutational_profile_counts(mutational_profile_counts, "IS", 'frobenius')
            if method == 'Z':
                _, _, contributing_signatures_A = decompose_mutational_profile_counts(mutational_profile_counts, "CLA", 'frobeniuszero')
                _, _, contributing_signatures_B = decompose_mutational_profile_counts(mutational_profile_counts, "CLB", 'frobeniuszero')
                _, _, contributing_signatures_cosmic = decompose_mutational_profile_counts(mutational_profile_counts, "IS", 'frobeniuszero')

            if method == 'D':
                _, _, contributing_signatures_cosmic = deconstruct_sigs(prefix + "/profiles/{}.profile".format(SAMPLE_ID), SAMPLE_ID)

            if method != 'D':
                oname = prefix + "/decompose/{}_{}-5.txt".format(SAMPLE_ID, method)
                with open(oname, 'w') as o:
                    for v in contributing_signatures_A:
                        o.write("{}\t{}\n".format(v['name'], v['score']))
                    # o.write("residuals\t{}\n".format(residuals_A))

                oname = prefix + "/decompose/{}_{}-10.txt".format(SAMPLE_ID, method)
                with open(oname, 'w') as o:
                    for v in contributing_signatures_B:
                        o.write("{}\t{}\n".format(v['name'], v['score']))
                    # o.write("residuals\t{}\n".format(residuals_B))

            oname = prefix + "/decompose/{}_{}-30.txt".format(SAMPLE_ID, method)
            with open(oname, 'w') as o:
                for v in contributing_signatures_cosmic:
                    o.write("{}\t{}\n".format(v['name'], v['score']))
                # o.write("residuals\t{}\n".format(residuals_cosmic))

    print(processing_stats)


if __name__ == '__main__':
    export_maf(
        prefix="/Users/gonceare/projects/Mutations/data/NCI60",
        maf_file="/Users/gonceare/projects/Mutations/data/NCI60/nci60.maf")

    # analyze_nci60_samples(prefix="/Users/gonceare/projects/Mutations/data/NCI60")

    # make_NCI60_report(prefix="/Users/gonceare/projects/Mutations/data/NCI60", n=30, decomposition='L')
    # make_NCI60_report(prefix="/Users/gonceare/projects/Mutations/data/NCI60", n=10, decomposition='L')
    # make_NCI60_report(prefix="/Users/gonceare/projects/Mutations/data/NCI60", n=5, decomposition='L')
