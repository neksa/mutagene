import csv
import logging
from collections import defaultdict, namedtuple

import twobitreader as tbr
from tqdm import tqdm

from mutagene.dna import complementary_nucleotide, nucleotides
from mutagene.io.context_stats import new_context_stats
from mutagene.io.maf_columns import normalize_header, report_malformed_rows, resolve_alleles
from mutagene.io.variant_filter import (
    DEFAULT_FILTER_COLUMN,
    check_filter_column,
    filter_value,
    passes_filter,
    report_filtered,
)

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

    methods = {"twobit": get_context_53_twobit}

    return methods[method](mutations, assembly)


def strip_line_terminator(line):
    """Remove the line terminator only, preserving empty trailing tab-delimited columns"""
    return line.rstrip("\r\n")


def read_auto_profile(muts, fmt, asm, keep_filtered=False, filter_column=DEFAULT_FILTER_COLUMN):
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

    # Sniffing leaves fmt as it found it, so an empty file never sets one and
    # concatenating it here was a TypeError before the format could be reported.
    logger.info(f"DATA FORMAT: {fmt}")

    if fmt not in ["MAF", "VCF"]:
        # An empty file, or one whose contents look like nothing we read. Report
        # an empty result rather than None: callers go on to read the counts out
        # of processing_stats, and None made an unreadable file a TypeError
        # deep in the caller instead of a message here.
        if fmt in (None, "", "AUTO"):
            logger.warning(
                "Could not tell whether this is a MAF or a VCF file. "
                "If the file is not empty, name the format with --input-format (-f)"
            )
        else:
            logger.warning(f"The dataformat [{fmt}] is not supported")
        return defaultdict(float), {"loaded": 0, "skipped": 0, "format": "unknown"}

    if fmt == "VCF":
        mutations, processing_stats = read_VCF_profile(mutations_lines, asm, keep_filtered)
    if fmt == "MAF":
        mutations, processing_stats = read_MAF_profile(
            mutations_lines, asm, keep_filtered, filter_column
        )

    return mutations, processing_stats


def read_MAF_profile(muts, asm, keep_filtered=False, filter_column=DEFAULT_FILTER_COLUMN):

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

    check_filter_column(header, filter_column, keep_filtered)

    N_loaded = N_skipped = 0

    context_stats = new_context_stats()
    raw_mutations = []
    skipped_rows = []
    n_filtered = 0
    for line_number, row in enumerate(reader, start=2):  # row 1 is the header
        # A row that cannot be read is skipped and counted rather than raising
        # out of the parser: one unusable row should not cost the whole file.
        try:
            data = MAF._make(row)
        except TypeError:
            skipped_rows.append((line_number, f"expected {len(header)} fields, got {len(row)}"))
            N_skipped += 1
            continue

        if not keep_filtered and not passes_filter(filter_value(data, filter_column)):
            n_filtered += 1
            continue

        chrom = getattr(data, "chromosome", None)
        if chrom is None:
            raise ValueError("Chromosome is not defined in MAF file")
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]

        try:
            pos = int(data.start_position)  # MAF POS START
            pos_end = int(data.end_position)  # MAF POS END
        except AttributeError as e:
            # A whole column is missing, which is a property of the file rather
            # than of this row, so there is no point reading any further.
            raise ValueError(
                f"Start_Position and End_Position are required in MAF file: {e}"
            ) from None
        except ValueError:
            skipped_rows.append((line_number, "position is not a number"))
            N_skipped += 1
            continue

        if pos != pos_end:
            continue

        x, y = resolve_alleles(data)  # MAF REF and variant allele
        if y is None:
            continue

        # skip if found unexpected nucleotide characters
        if len(set([x, y]) - set(nucleotides)) > 0:
            continue

        raw_mutations.append((chrom, pos, x, y))

    report_malformed_rows(skipped_rows, "MAF")
    report_filtered(n_filtered, "MAF")

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


def read_VCF_profile(muts, asm=None, keep_filtered=False):
    cn = complementary_nucleotide
    mutations = defaultdict(float)
    N_skipped = 0

    context_stats = new_context_stats()
    raw_mutations = []
    skipped_rows = []
    n_filtered = 0
    for i, line in enumerate(muts, start=1):
        if line.startswith("#"):
            continue
        if len(line) < 10:
            continue

        col_list = line.split()
        # ALT is the fifth column, so five are needed to read one. The old guard
        # asked for four and then indexed the fifth, so a four-column line raised
        # IndexError out of the parser instead of being skipped.
        if len(col_list) < 5:
            skipped_rows.append((i, f"expected at least 5 columns, got {len(col_list)}"))
            N_skipped += 1
            continue

        chrom = col_list[0]  # VCF CHROM
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]

        try:
            pos = int(col_list[1])  # VCF POS
        except ValueError:
            skipped_rows.append((i, f"position {col_list[1]!r} is not a number"))
            N_skipped += 1
            continue

        # FILTER is the seventh column. A file that stops before it has
        # recorded no filters, so there is nothing to honour.
        if not keep_filtered and len(col_list) > 6 and not passes_filter(col_list[6]):
            n_filtered += 1
            continue

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

    report_malformed_rows(skipped_rows, "VCF")
    report_filtered(n_filtered, "VCF")

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
