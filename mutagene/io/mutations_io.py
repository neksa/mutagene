from collections import defaultdict
from itertools import cycle

try:
    import twobitreader as tbr
except:
    print("twobitreader module required")

nucleotides = "ACGT"  # this order of nucleotides is important for reversing
complementary_nucleotide = dict(zip(nucleotides, reversed(nucleotides)))
complementary_nucleotide['N'] = 'N'

bases_dict = {"A": "A", "G": "G", "T": "T", "C": "C", "W": "AT", "S": "CG", "M": "AC", "K": "GT", "R": "AG", "Y": "CT",
              "B": "TCG", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ATGC"}

comp_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "W": "AT", "S": "CG", "K": "AC", "M": "GT", "Y": "AG",
             "R": "CT", "V": "TCG", "H": "AGT", "D": "ACT", "B": "ACG", "N": "ATGC"}

TWOBIT_GENOMES_PATH = '/net/pan1/mutagene/data/genomes/'


def get_context_twobit(mutations, assembly):
    """
    User twobitreader to get context of mutations
    """
    window_size = 50
    contexts = {}

    genomes_path = globals().get('TWOBIT_GENOMES_PATH', None)
    if not genomes_path:
        genomes_path = current_app.config['TWOBIT_GENOMES_PATH']
    if not genomes_path:
        return contexts

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

    if assembly not in twobit_files:
        return contexts
    else:
        twobit_file = genomes_path + "/" + twobit_files[assembly] + ".2bit"
        f = tbr.TwoBitFile(twobit_file)

        cn = complementary_nucleotide
        for (chrom, pos, x, y) in mutations:
            start = int(pos) - 1  # zero-based numbering
            chrom = str(chrom)

            chromosome = chrom if chrom.startswith('chr') else 'chr' + chrom
            chromosome = chromosome_name_mapping.get(chromosome, chromosome)

            if chromosome in f:
                try:
                    seq = f[chromosome][start - window_size:
                                        start + window_size + 1]
                    assert len(seq) == window_size * 2 + 1
                    seq = seq.upper()
                except Exception as e:
                    print("TwoBit exception", str(e), (chrom, pos, x, y))
            else:
                print("NO CHROM", chromosome)

            seq_with_coords = list(zip(
                cycle([chrom]),
                range(pos - window_size, pos + window_size + 1),
                seq,
                cycle("+")))

            nuc3 = seq_with_coords[window_size - 1][2]
            nuc = seq_with_coords[window_size][2]
            nuc5 = seq_with_coords[window_size + 1][2]

            if nuc != 'N' and nuc != x:
                if cn[nuc] == x:
                    nuc3 = cn[nuc5]
                    nuc5 = cn[nuc3]
                else:
                    print("{}:{}  {}>{}   {}[{}]{}".format(chromosome, pos, x, y, nuc5, nuc, nuc3))
                    nuc3 = nuc5 = 'N'
            contexts[(chrom, pos)] = (nuc5, nuc3), seq_with_coords
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
        'twobit': get_context_twobit
    }

    contexts = methods[method](mutations, assembly)
    return contexts


def read_mutations(muts, fmt=None, asm=None):
    mutations = None
    processing_stats = None
    if fmt is not None:
        fmt = fmt.upper()

    if fmt is None or fmt == "AUTO" or fmt == "":
        # for line in muts.split("\n"):
        for line in muts:
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
            continue

    print("DATA FORMAT:", fmt)

    if fmt is None or fmt == "AUTO" or fmt == "" or fmt == "TRI":
        try:
            mutations, mutations_with_context, processing_stats = read_trinucleotide(muts)
        except:
            pass

    if fmt == "VCF":
        mutations, mutations_with_context, processing_stats = read_VCF(muts, asm)
    if fmt == "MAF":
        mutations, mutations_with_context, processing_stats = read_MAF(muts, asm)
    if fmt == "PROFILE":
        mutations, mutations_with_context, processing_stats = read_PROFILE(muts)

    return mutations, mutations_with_context, processing_stats


def read_MAF(muts, asm=None):
    cn = complementary_nucleotide
    mutations = defaultdict(float)
    N_skipped = 0

    raw_mutations = []
    for i, line in enumerate(muts.split("\n")):
        if line.startswith("#"):
            continue
        if line.startswith("Hugo_Symbol"):
            continue
        if len(line) < 10:
            continue

        col_list = line.split()
        if len(col_list) < 13:
            continue

        try:
            # assembly_build = col_list[3]  # MAF ASSEMBLY
            # strand = col_list[7]    # MAF STRAND
            # ID = col_list[2]

            # chromosome is expected to be one or two number or one letter
            chrom = col_list[4]  # MAF CHROM
            # if chrom.lower().startswith("chr"):
            #     chrom = chrom[3:]

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

        if pos != pos_end:
            continue

        # skip if found unexpected nucleotide characters
        if len(set([x, y1, y2]) - set(nucleotides)) > 0:
            continue

        y = y1 if y1 != x else None
        y = y2 if y2 != x else y
        if y is None:
            continue

        raw_mutations.append((chrom, pos, x, y))

    mutations_with_context = []

    MAX = 100000000
    if len(raw_mutations) > 0:
        if len(raw_mutations) > MAX:
            raw_mutations = raw_mutations[:MAX]

        contexts = get_context_batch(raw_mutations, asm)

        if contexts is None:
            return None, None

        if len(contexts) == 0:
            return None, None

        for (chrom, pos, x, y) in raw_mutations:
            (p5, p3), seq_with_coords = contexts.get((chrom, pos), (("N", "N"), []))

            if len(set([p5, x, y, p3]) - set(nucleotides)) > 0:
                # print("Skipping invalid nucleotides")
                N_skipped += 1
                continue

            if x in "CT":
                mutations[p5 + p3 + x + y] += 1.0
            else:
                # complementary mutation
                mutations[cn[p3] + cn[p5] + cn[x] + cn[y]] += 1.0
            mutations_with_context.append((chrom, pos, x, y, seq_with_coords))

    N_loaded = int(sum(mutations.values()))
    processing_stats = {'loaded': N_loaded, 'skipped': N_skipped, 'format': 'MAF'}
    return mutations, mutations_with_context, processing_stats


def read_PROFILE(muts):
    mutations = defaultdict(float)

    for line in muts.split("\n"):
        if len(line) == 0:
            continue
        fields = line.upper().split()
        if len(fields) != 2:
            return None, None

        if len(fields[0]) != 7:
            return None, None

        if fields[0][1] != "[" or fields[0][3] != ">" or fields[0][5] != "]":
            return None, None

        p5, _, x, _, y, _, p3 = tuple(fields[0])

        if p5 not in nucleotides or p3 not in nucleotides or y not in nucleotides:
            return None, None

        if x not in "TC":
            return None, None

        try:
            f = float(fields[1])
        except:
            return None, None

        mutations[p5 + p3 + x + y] = f
        # print(p5, p3, x, y)

    processing_stats = {'loaded': 96, 'skipped': 0, 'format': 'mutational profile'}
    mutations_with_context = None
    return mutations, mutations_with_context, processing_stats


def read_trinucleotide(muts):
    N_skipped = 0
    cn = complementary_nucleotide
    mutations = defaultdict(float)

    for line in muts.split("\n"):
        if len(line.strip()) == 0:
            continue
        val = line.upper().strip().split(">")
        val = [v.strip() for v in val]
        if len(val) != 2:
            N_skipped += 1
            continue
        x3, y3 = val
        if len(x3) != 3 or len(y3) != 3:
            N_skipped += 1
            continue
        p5 = x3[0]
        x = x3[1]
        y = y3[1]
        p3 = x3[2]

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
    if N_loaded == 0:
        return None, None
    processing_stats = {'loaded': N_loaded, 'skipped': N_skipped, 'format': 'trinucleotides'}
    mutations_with_context = None
    return mutations, mutations_with_context, processing_stats


def read_VCF(muts, asm=None):
    cn = complementary_nucleotide
    mutations = defaultdict(float)
    N_skipped = 0
    N_skipped_indels = 0

    raw_mutations = []
    for line in muts.split("\n"):
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
    # print("RAW", raw_mutations)
    print("INDELS", N_skipped)

    MAX = 100000000
    mutations_with_context = []
    if len(raw_mutations) > 0:
        if len(raw_mutations) > MAX:
            raw_mutations = raw_mutations[:MAX]

        contexts = get_context_batch(raw_mutations, asm)
        # print("CONTEXTS", contexts)

        if contexts is None:
            return None, None

        if len(contexts) == 0:
            return None, None

        for (chrom, pos, x, y) in raw_mutations:
            (p5, p3), seq_with_coords = contexts.get((chrom, pos), (("N", "N"), []))
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

            mutations_with_context.append((chrom, pos, x, y, seq_with_coords))

    N_loaded = int(sum(mutations.values()))
    processing_stats = {'loaded': N_loaded, 'skipped': N_skipped, 'format': 'VCF'}
    return mutations, mutations_with_context, processing_stats
