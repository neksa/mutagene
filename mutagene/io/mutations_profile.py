import csv
from collections import namedtuple

import twobitreader as tbr

from tqdm import tqdm
from collections import defaultdict

from mutagene.dna import nucleotides, complementary_nucleotide

import logging
logger = logging.getLogger(__name__)


def get_context_53_twobit(mutations, twobit_file):
    """
    User twobitreader to get context of mutations
    """
    contexts = {}

    fname = twobit_file if twobit_file.endswith('.2bit') else twobit_file + '.2bit'
    f = tbr.TwoBitFile(fname)

    cn = complementary_nucleotide
    for (chrom, pos, x, y) in tqdm(mutations, leave=False):
        start = int(pos) - 1  # zero-based numbering
        chrom = str(chrom)
        chromosome = chrom if chrom.startswith('chr') else 'chr' + chrom

        nuc5 = 'N'
        nuc3 = 'N'
        nuc = 'N'
        if chromosome in f:
            try:
                seq = f[chromosome][start - 1: start + 2]  # +/- 1 nucleotide
                nuc5, nuc, nuc3 = tuple(seq.upper())
            except:
                nuc = 'N'
                # print(chromosome, x, nuc5, nuc, nuc3)
            if nuc != 'N' and nuc != x:
                if cn[nuc] == x:
                    nuc3 = cn[nuc5]
                    nuc5 = cn[nuc3]
                else:
                    nuc3 = nuc5 = 'N'
        else:
            # print("NO CHROM", chromosome)
            pass
        contexts[(chrom, pos)] = (nuc5, nuc3)
    return contexts


def get_context_batch(mutations, assembly, method='twobit'):
    """
        Get context for a list of mutations [(chrom, pos, x, y) ] format
    """
    if assembly is None:
        assembly = 38

    if method is None:
        method = 'twobit'

    methods = {
        # 'ensembl': get_context_ensembl,
        'twobit': get_context_53_twobit
    }

    contexts = methods[method](mutations, assembly)
    return contexts


def read_auto_profile(muts, fmt, asm):
    # if isinstance(muts, io.TextIOWrapper):
    #     print("OOOO")
    #     muts = muts.readlines()
    # print(muts)
    mutations = None
    processing_stats = None
    if fmt is not None:
        fmt = fmt.upper()

    if fmt is None or fmt == "AUTO" or fmt == 'auto' or fmt == "":
        mutations_lines = []
        for line in muts:
            mutations_lines.append(line.strip())
            if line.startswith("#version 2."):
                fmt = "MAF"
                break
            if len(line.strip()) == 0 or line.startswith("#"):
                continue
            tabs = line.split()
            if len(tabs) == 2 and "[" in tabs[0] and "]" in tabs[0]:
                fmt = "PROFILE"
                break
            if (len(tabs) == 3 and tabs[1] == ">") or (len(tabs) == 1 and ">" in tabs[0]):
                fmt = "TRI"
                break
            if len(tabs) > 3:
                fmt = "VCF"
                if tabs[0].lower().startswith("chr"):
                    break  # yes, it's VCF
                if tabs[4].lower().startswith("chr"):
                    fmt = "MAF"
                    break
        mutations_lines.extend(muts.readlines())
    else:
        mutations_lines = muts.readlines()

    logger.info("DATA FORMAT:" + fmt)

    if fmt not in ['MAF', 'VCF']:
        logger.warning("The dataformat [{}] is not supported".format(fmt))

    if fmt == "VCF":
        mutations, processing_stats = read_VCF_profile(mutations_lines, asm)
    if fmt == "MAF":
        mutations, processing_stats = read_MAF_profile(mutations_lines, asm)

    return mutations, processing_stats


def read_MAF_profile(muts, asm):

    cn = complementary_nucleotide
    mutations = defaultdict(float)
    N_skipped = 0

    try:
        reader = csv.reader((row for row in muts if not row.startswith('#')), delimiter='\t')
        # get names from column headers
        header = next(reader)
        # MAF = namedtuple("MAF", map(str.lower, next(reader)))
        header = tuple(map(lambda s: s.lower().replace('.', '_'), header))
        # print(header)
        MAF = namedtuple("MAF", header, rename=True)
    except ValueError:
        # raise
        logger.warning("MAF format not recognized")
        return mutations, {}

    N_loaded = N_skipped = 0

    raw_mutations = []
    for data in map(MAF._make, reader):
        # assembly_build = col_list[3]  # MAF ASSEMBLY
        # strand = col_list[7]    # MAF STRAND
        # ID = col_list[2]

        # chromosome is expected to be one or two number or one letter
        chrom = data.chromosome  # MAF CHROM
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]
        # if len(chrom) == 2 and chrom[1] not in "0123456789":
        #     chrom = chrom[0]

        pos = int(data.start_position)    # MAF POS START
        pos_end = int(data.end_position)  # MAF POS END
        x = data.reference_allele         # MAF REF

        y1 = data.tumor_seq_allele1       # MAF ALT1
        y2 = data.tumor_seq_allele2       # MAF ALT2

        if pos != pos_end:
            continue

        # skip if found unexpected nucleotide characters
        if len(set([x, y1, y2]) - set(nucleotides)) > 0:
            continue
        # print("OK")

        y = y1 if y1 != x else None
        y = y2 if y2 != x else y
        if y is None:
            continue

        raw_mutations.append((chrom, pos, x, y))

    # MAX = 100000000
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
                N_skipped += 1
                continue

            if x in "CT":
                mutations[p5 + p3 + x + y] += 1.0
            else:
                # complementary mutation
                mutations[cn[p3] + cn[p5] + cn[x] + cn[y]] += 1.0

    N_loaded = int(sum(mutations.values()))
    processing_stats = {'loaded': N_loaded, 'skipped': N_skipped, 'format': 'MAF'}
    return mutations, processing_stats


def read_VCF_profile(muts, asm=None):
    cn = complementary_nucleotide
    mutations = defaultdict(float)
    N_skipped = 0

    raw_mutations = []
    for i, line in enumerate(muts):
        if line.startswith("#"):
            continue
        if len(line) < 10:
            continue

        col_list = line.split()
        if len(col_list) < 4:
            continue

        # ID = col_list[2]
        # chromosome is expected to be one or two number or one letter
        chrom = col_list[0]  # VCF CHROM
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]
        # if len(chrom) == 2 and chrom[1] not in "0123456789":
        #     chrom = chrom[0]

        pos = int(col_list[1])  # VCF POS
        x = col_list[3]         # VCF REF
        y = col_list[4]         # VCF ALT

        # if multiple REF or ALT alleles are given, ignore mutation entry (could mean seq error, could mean deletion)
        if len(x) != 1:
            N_skipped += 1
            continue
        if len(y) != 1:
            N_skipped += 1
            continue

        raw_mutations.append((chrom, pos, x, y))

    # MAX = 1e7
    if len(raw_mutations) > 0:
        # if len(raw_mutations) > MAX:
        #     raw_mutations = raw_mutations[:MAX]

        contexts = get_context_batch(raw_mutations, asm)
        # print("CONTEXTS", contexts)

        if contexts is None:
            return None, None

        if len(contexts) == 0:
            return None, None

        for (chrom, pos, x, y) in raw_mutations:
            p5, p3 = contexts.get((chrom, pos), ("N", "N"))
            # print("RESULT: {} {}".format(p5, p3))

            if len(set([p5, x, y, p3]) - set(nucleotides)) > 0:
                # print(chrom, pos, p5, p3, x)
                # print("Skipping invalid nucleotides")
                N_skipped += 1
                continue

            if x in "CT":
                mutations[p5 + p3 + x + y] += 1.0
            else:
                # complementary mutation
                mutations[cn[p3] + cn[p5] + cn[x] + cn[y]] += 1.0

    N_loaded = int(sum(mutations.values()))
    processing_stats = {'loaded': N_loaded, 'skipped': N_skipped, 'format': 'VCF'}
    return mutations, processing_stats
