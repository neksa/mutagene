import csv
from collections import namedtuple
import twobitreader as tbr
from tqdm import tqdm
from collections import defaultdict
from itertools import cycle

from mutagene.dna import chromosome_name_mapping
from mutagene.motifs.motifs import nucleotides, complementary_nucleotide

import logging
logger = logging.getLogger(__name__)


def get_context_twobit_window(mutations, twobit_file, window_size=50):
    """
    User twobitreader to get context of mutations
    """
    contexts = {}

    fname = twobit_file if twobit_file.endswith('.2bit') else twobit_file + '.2bit'
    f = tbr.TwoBitFile(fname)

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
                logger.warning("TwoBit exception: ", str(e), (chrom, pos, x, y))
        else:
            logger.warning("NO CHROM", chromosome)

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


def read_MAF_with_context_window(infile, asm=None):
    cn = complementary_nucleotide
    mutations = defaultdict(lambda: defaultdict(float))
    N_skipped = 0

    processing_stats = {'loaded': 0, 'skipped': 0, 'nsamples': 0, 'format': 'unknown'}
    if not infile:
        logger.warning("No input file")
        return mutations, {}, processing_stats

    try:
        reader = csv.reader((row for row in infile if not row.startswith('#')), delimiter='\t')
        # get names from column headers
        header = next(reader)
        header = tuple(map(lambda s: s.lower().replace('.', '_'), header))
        # print(header)
        MAF = namedtuple("MAF", header, rename=True)
    except ValueError:
        raise
        logger.warning("MAF format not recognized")
        return mutations, {}, processing_stats

    raw_mutations = defaultdict(list)
    # for line in tqdm(infile):
    for data in tqdm(map(MAF._make, reader)):
        try:
            # assembly_build = col_list[3]  # MAF ASSEMBLY
            # strand = col_list[7]    # MAF STRAND
            # ID = col_list[2]

            # chromosome is expected to be one or two number or one letter
            chrom = data.chromosome  # MAF CHROM
            # if chrom.lower().startswith("chr"):
            #     chrom = chrom[3:]

            # if len(chrom) == 2 and chrom[1] not in "0123456789":
            #     chrom = chrom[0]

            pos = int(data.start_position)
            # pos_end = int(col_list[6])  # MAF POS END
            x = data.reference_allele            # MAF REF
            y1 = data.tumor_seq_allele1           # MAF ALT1
            y2 = data.tumor_seq_allele2           # MAF ALT2
            sample = data.tumor_sample_barcode       # Tumor barcode
        except:
            # raise
            continue

        # if pos != pos_end:
        #     continue

        # skip if found unexpected nucleotide characters
        if len(set([x, y1, y2]) - set(nucleotides)) > 0:
            continue

        y = y1 if y1 != x else None
        y = y2 if y2 != x else y
        if y is None:
            continue

        raw_mutations[sample].append((chrom, pos, x, y))

    mutations_with_context = defaultdict(list)

    for sample, sample_mutations in raw_mutations.items():
        if len(sample_mutations) > 0:
            contexts = get_context_twobit_window(sample_mutations, asm)

            if contexts is None or len(contexts) == 0:
                return None, None

            for (chrom, pos, x, y) in sample_mutations:
                (p5, p3), seq_with_coords = contexts.get((chrom, pos), (("N", "N"), []))

                if len(set([p5, x, y, p3]) - set(nucleotides)) > 0:
                    # print("Skipping invalid nucleotides")
                    N_skipped += 1
                    continue

                if x in "CT":
                    mutations[sample][p5 + p3 + x + y] += 1.0
                else:
                    # complementary mutation
                    mutations[sample][cn[p3] + cn[p5] + cn[x] + cn[y]] += 1.0
                mutations_with_context[sample].append((chrom, pos, x, y, seq_with_coords))

    N_loaded = 0
    for sample, sample_mutations in mutations.items():
        N_loaded += int(sum(sample_mutations.values()))
    processing_stats = {
        'loaded': N_loaded,
        'skipped': N_skipped,
        'nsamples': len(mutations.keys()),
        'format': 'MAF'
    }
    return mutations, mutations_with_context, processing_stats


def read_VCF_with_context_window(muts, asm=None):
    cn = complementary_nucleotide
    mutations = defaultdict(float)
    N_skipped = 0
    # N_skipped_indels = 0

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
    # print("INDELS", N_skipped)

    mutations_with_context = []
    if len(raw_mutations) > 0:
        contexts = get_context_twobit_window(raw_mutations, asm)
        if contexts is None or len(contexts) == 0:
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
