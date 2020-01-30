import csv
from collections import namedtuple
import twobitreader as tbr
from tqdm import tqdm
from collections import defaultdict
from itertools import cycle

from mutagene.dna import chromosome_name_mapping
from mutagene.motifs import nucleotides, complementary_nucleotide

import logging
logger = logging.getLogger(__name__)


def get_context_twobit_window(mutations, twobit_file, window_size):
    """
    User twobitreader to get context of mutations
    """
    if window_size is None:
        window_size = 50

    contexts = {}

    fname = twobit_file if twobit_file.endswith('.2bit') else twobit_file + '.2bit'
    f = tbr.TwoBitFile(fname)

    cn = complementary_nucleotide
    for (chrom, pos, transcript_strand, x, y) in mutations:
        start = int(pos) - 1  # 2bit uses zero-based numbering
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
                logger.warning("TwoBit exception while reading the genome in {}:{}: {}".format(chrom, pos, e))
                continue
        else:
            logger.warning("Chromosome {} not found in 2bit file. Consider renaming it or using a different genome assembly".format(chromosome))
            continue

        strand = transcript_strand
        seq_with_coords = list(zip(
            cycle([chrom]),
            range(pos - window_size, pos + window_size + 1),
            seq,
            cycle([strand])))

        assert len(seq_with_coords) == len(seq)

        nuc5 = seq_with_coords[window_size - 1][2]
        nuc = seq_with_coords[window_size][2]
        nuc3 = seq_with_coords[window_size + 1][2]

        if nuc != 'N' and nuc != x:
            if cn[nuc] == x:
                nuc3 = cn[nuc5]
                nuc5 = cn[nuc3]
                # print('debug: complementary REF sequence detected')
            else:
                # print("{}:{}  {}>{}   {}[{}]{}".format(chromosome, pos, x, y, nuc5, nuc, nuc3))
                nuc3 = nuc5 = 'N'
            logger.warning(
                "REF allele does not match the genomic sequence in {}:{} {}!={}. Multiple errors could mean wrong genome assembly choice".format(
                    chromosome, pos, x, nuc))
        contexts[(chrom, pos)] = (nuc5, nuc3), seq_with_coords
    return contexts


def read_MAF_with_context_window(infile, asm, window_size):
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
        logger.warning("MAF format not recognized")
        raise
        return mutations, {}, processing_stats

    raw_mutations = defaultdict(list)
    # for line in tqdm(infile):
    for data in tqdm(map(MAF._make, reader), leave=False):
        # assembly_build = col_list[3]  # MAF ASSEMBLY

        # chromosome is expected to be one or two number or one letter
        if hasattr(data, 'chromosome'):
            chrom = data.chromosome  # MAF CHROM
        else:
            raise ValueError('Chromosome is not defined in MAF file')

        if hasattr(data, 'tumor_sample_barcode'):
            sample = data.tumor_sample_barcode       # Tumor barcode
        elif hasattr(data, 'sample_id'):
            sample = data.sample_id
        else:
            raise ValueError("Sample ID is not defined in MAF file")

        if hasattr(data, 'reference_allele'):
            x = data.reference_allele            # MAF REF
        else:
            raise ValueError('Reference allele is not defined in MAF file')

        if hasattr(data, 'variant_allele'):
            y = data.variant_allele
        elif hasattr(data, 'tumor_seq_allele1') and hasattr(data, 'tumor_seq_allele2'):
            y1 = data.tumor_seq_allele1           # MAF ALT1
            y2 = data.tumor_seq_allele2           # MAF ALT2
            y = y1 if y1 != x else None
            y = y2 if y2 != x else y
        else:
            raise ValueError('Variant allele is not defined in MAF file')

        if y is None:
            continue
        # skip if found unexpected nucleotide characters
        if len(set([x, y]) - set(nucleotides)) > 0:
            continue

        if hasattr(data, 'start_position'):
            try:
                pos = int(data.start_position)
            except ValueError:
                raise ValueError('Start position is not a number in MAF file')
        else:
            raise ValueError('Start position is not defined in MAF file')

        if hasattr(data, 'transcript_strand'):
            transcript_strand = data.transcript_strand
            # GDC uses 1 and -1
            if transcript_strand == '+':
                pass
            if transcript_strand == '-':
                pass
            elif transcript_strand == '1':
                transcript_strand = '+'
            elif transcript_strand == '-1':
                transcript_strand = '-'
            elif transcript_strand == '':
                transcript_strand = '+'  # default value
            else:
                raise ValueError('Unexpected value of transcription_strand in MAF file')
        else:
            # this is an incorrect assumption about transcription strand
            transcript_strand = '+'

        raw_mutations[sample].append((chrom, pos, transcript_strand, x, y))

    mutations_with_context = defaultdict(list)

    for sample, sample_mutations in raw_mutations.items():
        if len(sample_mutations) > 0:
            contexts = get_context_twobit_window(sample_mutations, asm, window_size)

            if contexts is None or len(contexts) == 0:
                return None, None

            for (chrom, pos, transcript_strand, x, y) in sample_mutations:
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
                mutations_with_context[sample].append((chrom, pos, transcript_strand, x, y, seq_with_coords))

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


def read_VCF_with_context_window(infile, asm, window_size):
    cn = complementary_nucleotide
    mutations = defaultdict(lambda: defaultdict(float))
    mutations_with_context = defaultdict(list)
    raw_mutations = defaultdict(list)

    N_skipped = 0
    # N_skipped_indels = 0

    sample = 'VCF'

    for line in infile:
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

        transcript_strand = '+'
        raw_mutations[sample].append((chrom, pos, transcript_strand, x, y))
    # print("RAW", raw_mutations)
    # print("INDELS", N_skipped)

    for sample, sample_mutations in raw_mutations.items():
        if len(sample_mutations) > 0:
            contexts = get_context_twobit_window(sample_mutations, asm, window_size)
            if contexts is None or len(contexts) == 0:
                return None, None

            for (chrom, pos, transcript_strand, x, y) in sample_mutations:
                (p5, p3), seq_with_coords = contexts.get((chrom, pos), (("N", "N"), []))
                # print("RESULT: {} {}".format(p5, p3))

                if len(set([p5, x, y, p3]) - set(nucleotides)) > 0:
                    # print(chrom, pos, p5, p3, x)
                    # print("Skipping invalid nucleotides")
                    N_skipped += 1
                    continue

                if x in "CT":
                    mutations[sample][p5 + p3 + x + y] += 1.0
                else:
                    # complementary mutation
                    mutations[sample][cn[p3] + cn[p5] + cn[x] + cn[y]] += 1.0

                mutations_with_context[sample].append((chrom, pos, transcript_strand, x, y, seq_with_coords))

    N_loaded = 0
    for sample, sample_mutations in mutations.items():
        N_loaded += int(sum(sample_mutations.values()))

    nsamples = len(mutations.keys())
    processing_stats = {
        'loaded': N_loaded,
        'skipped': N_skipped,
        'format': 'VCF',
        'nsamples': nsamples
    }
    return mutations, mutations_with_context, processing_stats
