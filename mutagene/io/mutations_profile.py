import csv
import logging
from collections import defaultdict, namedtuple

import twobitreader as tbr
from tqdm import tqdm

from mutagene.dna import complementary_nucleotide, nucleotides
from mutagene.io.context_stats import new_context_stats
from mutagene.io.maf_columns import normalize_header, resolve_alleles

logger = logging.getLogger(__name__)


def get_context_53_twobit(mutations, twobit_file):
    """
    User twobitreader to get context of mutations

    returns contexts, context_stats
    """
    contexts = {}
    stats = new_context_stats()

    fname = twobit_file if twobit_file.endswith(".2bit") else twobit_file + ".2bit"
    f = tbr.TwoBitFile(fname)

    cn = complementary_nucleotide
    for chrom, pos, x, y in tqdm(mutations, leave=False):
        start = int(pos) - 1  # zero-based numbering
        chrom = str(chrom)
        chromosome = chrom if chrom.startswith("chr") else "chr" + chrom

        nuc5 = "N"
        nuc3 = "N"
        nuc = "N"
        if chromosome in f:
            try:
                seq = f[chromosome][start - 1 : start + 2]  # +/- 1 nucleotide
                nuc5, nuc, nuc3 = tuple(seq.upper())
            except (ValueError, KeyError, IndexError):
                nuc = "N"
                stats["read_errors"] += 1
            if nuc != "N" and nuc != x:
                if cn[nuc] == x:
                    # REF is reported on the strand opposite the assembly, so the
                    # flanking bases swap sides as well as complement. Assigning
                    # them one at a time made the second read the value the first
                    # had just written (#101).
                    nuc5, nuc3 = cn[nuc3], cn[nuc5]
                    stats["reverse_strand_ref"] += 1
                else:
                    # REF matches neither strand: the wrong-assembly signal. This
                    # path used to drop the mutation silently, which is half the
                    # reason the webapp's mismatch warning never fired (#99).
                    nuc3 = nuc5 = "N"
                    stats["ref_mismatches"] += 1
        else:
            stats["chromosome_not_found"] += 1
        contexts[(chrom, pos)] = (nuc5, nuc3)
    return contexts, stats


def get_context_batch(mutations, assembly, method="twobit"):
    """
    Get context for a list of mutations [(chrom, pos, x, y) ] format

    returns contexts, context_stats
    """
    if assembly is None:
        assembly = 38

    if method is None:
        method = "twobit"

    methods = {
        # 'ensembl': get_context_ensembl,
        "twobit": get_context_53_twobit
    }

    return methods[method](mutations, assembly)


def strip_line_terminator(line):
    """Remove the line terminator only, preserving empty trailing tab-delimited columns"""
    return line.rstrip("\r\n")


def read_auto_profile(muts, fmt, asm):
    mutations = None
    processing_stats = None
    if fmt is not None:
        fmt = fmt.upper()

    if fmt is None or fmt == "AUTO" or fmt == "auto" or fmt == "":
        mutations_lines = []
        for line in muts:
            mutations_lines.append(strip_line_terminator(line))
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
                # Check for MAF-specific column names
                line_lower = line.lower()
                maf_indicators = [
                    "hugo_symbol",
                    "start_position",
                    "end_position",
                    "reference_allele",
                    "tumor_seq_allele",
                ]
                if any(indicator in line_lower for indicator in maf_indicators):
                    fmt = "MAF"
                    break
                # Fall back to position-based detection
                fmt = "VCF"
                if tabs[0].lower().startswith("chr"):
                    break  # yes, it's VCF
                if len(tabs) > 4 and tabs[4].lower().startswith("chr"):
                    fmt = "MAF"
                    break
        mutations_lines.extend(strip_line_terminator(line) for line in muts.readlines())
    else:
        mutations_lines = [strip_line_terminator(line) for line in muts.readlines()]

    logger.info("DATA FORMAT:" + fmt)

    if fmt not in ["MAF", "VCF"]:
        logger.warning(f"The dataformat [{fmt}] is not supported")

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
        reader = csv.reader((row for row in muts if not row.startswith("#")), delimiter="\t")
        header = next(reader)
        header = normalize_header(header)
        MAF = namedtuple("MAF", header, rename=True)
    except ValueError:
        # raise
        logger.warning("MAF format not recognized")
        return mutations, {}

    N_loaded = N_skipped = 0

    context_stats = new_context_stats()
    raw_mutations = []
    for data in map(MAF._make, reader):
        chrom = data.chromosome  # MAF CHROM
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]

        pos = int(data.start_position)  # MAF POS START
        pos_end = int(data.end_position)  # MAF POS END

        if pos != pos_end:
            continue

        x, y = resolve_alleles(data)  # MAF REF and variant allele
        if y is None:
            continue

        # skip if found unexpected nucleotide characters
        if len(set([x, y]) - set(nucleotides)) > 0:
            continue

        raw_mutations.append((chrom, pos, x, y))

    if len(raw_mutations) > 0:
        contexts, context_stats = get_context_batch(raw_mutations, asm)

        if contexts is None:
            return None, None

        if len(contexts) == 0:
            return None, None

        for chrom, pos, x, y in raw_mutations:
            p5, p3 = contexts.get((chrom, pos), ("N", "N"))

            if len(set([p5, x, y, p3]) - set(nucleotides)) > 0:
                N_skipped += 1
                continue

            if x in "CT":
                mutations[p5 + p3 + x + y] += 1.0
            else:
                # complementary mutation
                mutations[cn[p3] + cn[p5] + cn[x] + cn[y]] += 1.0

    N_loaded = int(sum(mutations.values()))
    processing_stats = {
        "loaded": N_loaded,
        "skipped": N_skipped,
        "format": "MAF",
        **context_stats,
    }
    return mutations, processing_stats


def read_VCF_profile(muts, asm=None):
    cn = complementary_nucleotide
    mutations = defaultdict(float)
    N_skipped = 0

    context_stats = new_context_stats()
    raw_mutations = []
    for i, line in enumerate(muts):
        if line.startswith("#"):
            continue
        if len(line) < 10:
            continue

        col_list = line.split()
        if len(col_list) < 4:
            continue

        chrom = col_list[0]  # VCF CHROM
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]

        pos = int(col_list[1])  # VCF POS
        x = col_list[3]  # VCF REF
        y = col_list[4]  # VCF ALT

        # if multiple REF or ALT alleles are given, ignore mutation entry (could mean seq error, could mean deletion)
        if len(x) != 1:
            N_skipped += 1
            continue
        if len(y) != 1:
            N_skipped += 1
            continue

        raw_mutations.append((chrom, pos, x, y))

    if len(raw_mutations) > 0:
        contexts, context_stats = get_context_batch(raw_mutations, asm)

        if contexts is None:
            return None, None

        if len(contexts) == 0:
            return None, None

        for chrom, pos, x, y in raw_mutations:
            p5, p3 = contexts.get((chrom, pos), ("N", "N"))

            if len(set([p5, x, y, p3]) - set(nucleotides)) > 0:
                N_skipped += 1
                continue

            if x in "CT":
                mutations[p5 + p3 + x + y] += 1.0
            else:
                # complementary mutation
                mutations[cn[p3] + cn[p5] + cn[x] + cn[y]] += 1.0

    N_loaded = int(sum(mutations.values()))
    processing_stats = {
        "loaded": N_loaded,
        "skipped": N_skipped,
        "format": "VCF",
        **context_stats,
    }
    return mutations, processing_stats
