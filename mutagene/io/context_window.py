import csv
import logging
from collections import defaultdict, namedtuple
from itertools import cycle

import twobitreader as tbr
from tqdm import tqdm

from mutagene.dna import chromosome_name_mapping
from mutagene.io.context_stats import merge_context_stats, new_context_stats
from mutagene.io.maf_columns import (
    normalize_header,
    report_malformed_rows,
    resolve_alleles,
    resolve_sample,
)
from mutagene.io.variant_filter import passes_filter, report_filtered
from mutagene.motifs import complementary_nucleotide, nucleotides

logger = logging.getLogger(__name__)


def get_context_twobit_window(mutations, twobit_file, window_size):
    """
    User twobitreader to get context of mutations
    It's a special data structure.
    contexts[(chrom, pos)] = (nuc5, nuc3), seq_with_coords
    where seq_with_coords = [(chrom, pos, nucleotide, strand)]

    returns contexts, context_stats
    """
    if window_size is None:
        window_size = 50

    contexts = {}
    stats = new_context_stats()

    fname = twobit_file if twobit_file.endswith(".2bit") else twobit_file + ".2bit"

    try:
        f = tbr.TwoBitFile(fname)
    except FileNotFoundError:
        raise FileNotFoundError(f"The 2bit genome assembly file, {fname}, was not found!")
    except Exception:
        raise RuntimeError(
            f"An error occurred while reading the 2bit genome assembly file, {fname}!"
        )

    cn = complementary_nucleotide
    for chrom, pos, transcript_strand, x, y in mutations:
        start = int(pos) - 1  # 2bit uses zero-based numbering
        chrom = str(chrom)

        chromosome = chrom if chrom.startswith("chr") else "chr" + chrom
        chromosome = chromosome_name_mapping.get(chromosome, chromosome)

        if chromosome in f:
            try:
                seq = f[chromosome][start - window_size : start + window_size + 1]
                assert len(seq) == window_size * 2 + 1
                seq = seq.upper()
            except Exception as e:
                logger.warning(f"TwoBit exception while reading the genome in {chrom}:{pos}: {e}")
                stats["read_errors"] += 1
                continue
        else:
            logger.warning(
                f"Chromosome {chromosome} not found in 2bit file. Consider renaming it or using a different genome assembly"
            )
            stats["chromosome_not_found"] += 1
            continue

        strand = transcript_strand
        seq_with_coords = list(
            zip(
                cycle([chrom]),
                range(pos - window_size, pos + window_size + 1),
                seq,
                cycle([strand]),
            )
        )

        assert len(seq_with_coords) == len(seq)

        nuc5 = seq_with_coords[window_size - 1][2]
        nuc = seq_with_coords[window_size][2]
        nuc3 = seq_with_coords[window_size + 1][2]

        if nuc != "N" and nuc != x:
            if cn[nuc] == x:
                # REF reported on the strand opposite the assembly. Normal, and
                # not evidence of a wrong genome, so counted separately.
                nuc5, nuc3 = cn[nuc3], cn[nuc5]
                stats["reverse_strand_ref"] += 1
            else:
                # REF matches neither strand: this is the wrong-assembly signal.
                # print("{}:{}  {}>{}   {}[{}]{}".format(chromosome, pos, x, y, nuc5, nuc, nuc3))
                nuc3 = nuc5 = "N"
                stats["ref_mismatches"] += 1
                logger.warning(
                    f"REF allele does not match the genomic sequence in {chromosome}:{pos} {x}!={nuc}. Multiple errors could mean wrong genome assembly choice"
                )
        contexts[(chrom, pos)] = (nuc5, nuc3), seq_with_coords
    return contexts, stats


def _assemble_mutations(raw_mutations, asm, window_size):
    """Attach genomic context to raw mutations, one sample at a time.

    Returns ``(mutations, mutations_with_context, n_skipped, context_stats)``.

    A sample that yields no context at all -- because its chromosome names are
    absent from the chosen assembly, say -- is dropped and the remaining samples
    are still processed. Returning from the whole function there (#100) meant a
    single bad sample discarded every sample after it and the caller saw
    ``loaded: 0``, as though the file were empty, with sample order deciding how
    much was lost. A genome that genuinely cannot be read still raises from
    ``get_context_twobit_window``.
    """
    cn = complementary_nucleotide
    mutations = defaultdict(lambda: defaultdict(float))
    mutations_with_context = defaultdict(list)
    samples_without_context = []
    context_stats = new_context_stats()
    n_skipped = 0

    for sample, sample_mutations in raw_mutations.items():
        if len(sample_mutations) == 0:
            continue

        contexts, sample_stats = get_context_twobit_window(sample_mutations, asm, window_size)
        merge_context_stats(context_stats, sample_stats)

        if contexts is None or len(contexts) == 0:
            samples_without_context.append(sample)
            n_skipped += len(sample_mutations)
            continue

        for chrom, pos, transcript_strand, x, y in sample_mutations:
            (p5, p3), seq_with_coords = contexts.get((chrom, pos), (("N", "N"), []))

            if len(set([p5, x, y, p3]) - set(nucleotides)) > 0:
                n_skipped += 1
                continue

            if x in "CT":
                mutations[sample][p5 + p3 + x + y] += 1.0
            else:
                # complementary mutation
                mutations[sample][cn[p3] + cn[p5] + cn[x] + cn[y]] += 1.0
            mutations_with_context[sample].append(
                (chrom, pos, transcript_strand, x, y, seq_with_coords)
            )

    if samples_without_context:
        shown = ", ".join(map(str, samples_without_context[:5]))
        if len(samples_without_context) > 5:
            shown += f" and {len(samples_without_context) - 5} more"
        logger.warning(
            f"No genomic context found for {len(samples_without_context)} sample(s): {shown}. "
            "Check that the chromosome names match the chosen genome assembly. "
            "The remaining samples were processed"
        )

    return mutations, mutations_with_context, n_skipped, context_stats


def read_TCGI_with_context_window(infile, asm, window_size):
    """
    Tabular file; no particular column order required but must contain header line with four mandatory column names:
    (CHR, POS, REF, ALT) corresponding to the chromosome, position, reference and alternate allele columns, respectively
    • Optionally, if column name SAMPLE exists, the column is used as sample names, otherwise it is assumed that all variants are in the same sample
    • Optionally, if column name STRAND exists, the column is used as strandedness of variants
    (possible values: 1 or + for forward strand, -1 or - for reverse strand), otherwise it is assumed that all variants are described on the forward strand
    • The file can contain comment lines starting with #
    • In addition to the CHR, POS, REF and ALT columns, can have extra columns that will also appear in the final TCGI
    output table
    • Coordinates must follow the GRCh37 genome build, chromosome names are accepted either with or without the ‘chr’ prefix
    • Alleles can follow both the Ensembl or VCF convention (e.g. for a deletion, both ‘ATCA to A’ or ‘TCA to –’ forms are accepted)

    https://www.cancergenomeinterpreter.org/images/inputformats.pdf

    returns mutations, mutations_with_context, processing_stats
    """
    N_skipped = 0

    processing_stats = {"loaded": 0, "skipped": 0, "nsamples": 0, "format": "unknown"}
    if not infile:
        logger.warning("No input file")
        return defaultdict(lambda: defaultdict(float)), {}, processing_stats

    try:
        reader = csv.reader((row for row in infile if not row.startswith("#")), delimiter="\t")
        # get names from column headers
        header = next(reader)
        header = normalize_header(header)
        # print(header)
        TCGI = namedtuple("TCGI", header, rename=True)
    except ValueError:
        logger.warning("TCGI format not recognized")
        raise

    raw_mutations = defaultdict(list)
    # for line in tqdm(infile):
    for data in tqdm(map(TCGI._make, reader), leave=False):
        # chromosome is expected to be one or two number or one letter.
        # The documented column name is CHR; CHROM is accepted because the code
        # only ever looked for that, so files written to the documentation were
        # rejected and files written to the code were undocumented.
        chrom = getattr(data, "chr", None)
        if chrom is None:
            chrom = getattr(data, "chrom", None)
        if chrom is None:
            raise ValueError("Chromosome (CHR) is not defined in TCGI file")

        if hasattr(data, "sample"):
            sample = data.sample
        else:
            sample = "TCGI"

        if hasattr(data, "ref"):
            x = data.ref
        else:
            raise ValueError("Reference allele (REF) is not defined in TCGI file")

        if hasattr(data, "alt"):
            y = data.alt
        else:
            raise ValueError("Variant allele is not defined in TCGI file")

        if y is None:
            continue
        # skip if found unexpected nucleotide characters
        if len(set([x, y]) - set(nucleotides)) > 0:
            continue

        if hasattr(data, "pos"):
            try:
                pos = int(data.pos)
            except ValueError:
                raise ValueError("Start position is not a number in TCGI file")
        else:
            raise ValueError("Start position is not defined in TCGI file")

        transcript_strand = "+"
        raw_mutations[sample].append((chrom, pos, transcript_strand, x, y))

    mutations, mutations_with_context, n_skipped, context_stats = _assemble_mutations(
        raw_mutations, asm, window_size
    )
    N_skipped += n_skipped

    N_loaded = 0
    for sample, sample_mutations in mutations.items():
        N_loaded += int(sum(sample_mutations.values()))

    processing_stats = {
        "loaded": N_loaded,
        "skipped": N_skipped,
        "nsamples": len(mutations.keys()),
        "format": "TCGI",
        **context_stats,
    }
    return mutations, mutations_with_context, processing_stats


_SUPPORTED_FORMATS = {"MAF", "VCF", "TCGI"}


class _LineNumberTracker:
    """
    Iterates over input lines dropping comments, remembering the source line number
    of the line handed to csv.reader last, so that malformed rows can be reported.
    """

    def __init__(self, lines):
        self.lines = lines
        self.line_number = 0

    def __iter__(self):
        for line_number, line in enumerate(self.lines, start=1):
            self.line_number = line_number
            if line.startswith("#"):
                continue
            yield line


# GDC uses 1 and -1, an empty value means the strand is unknown
_TRANSCRIPT_STRAND_VALUES = {"+": "+", "-": "-", "1": "+", "-1": "-", "": "+"}


def _normalize_transcript_strand(value):
    """Map a MAF transcript_strand value to '+' or '-', raising ValueError on unknown values"""
    try:
        return _TRANSCRIPT_STRAND_VALUES[value]
    except KeyError:
        raise ValueError(f"unexpected value of transcript_strand {value!r}") from None


def read_mutations(file_format, *args, **kwargs):
    """Wrapper for read_X_with_context_window"""
    if file_format not in _SUPPORTED_FORMATS:
        raise ValueError(
            f"Unsupported file format: {file_format}. Valid options: {', '.join(sorted(_SUPPORTED_FORMATS))}"
        )
    function_name = f"read_{file_format}_with_context_window"
    return globals()[function_name](*args, **kwargs)


def read_MAF_with_context_window(infile, asm, window_size, keep_filtered=False):
    """
    Read MAF file and extract context of mutations for assembly asm and window +/- window_size around each mutation
    MAF format description: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/

    returns mutations, mutations_with_context, processing_stats
    """
    N_skipped = 0

    processing_stats = {"loaded": 0, "skipped": 0, "nsamples": 0, "format": "unknown"}
    if not infile:
        logger.warning("No input file")
        return defaultdict(lambda: defaultdict(float)), {}, processing_stats

    tracker = _LineNumberTracker(infile)
    try:
        reader = csv.reader(tracker, delimiter="\t")
        # get names from column headers
        header = next(reader)
        header = normalize_header(header)
        # print(header)
        MAF = namedtuple("MAF", header, rename=True)
    except ValueError:
        logger.warning("MAF format not recognized")
        raise

    raw_mutations = defaultdict(list)
    skipped_rows = []
    n_filtered = 0
    # for line in tqdm(infile):
    for row in tqdm(reader, leave=False):
        if not any(field.strip() for field in row):
            # empty line, carries no mutation
            continue

        try:
            data = MAF._make(row)
        except TypeError:
            # wrong number of fields in this row: skip it and keep reading the file
            skipped_rows.append(
                (tracker.line_number, f"expected {len(header)} fields, got {len(row)}")
            )
            N_skipped += 1
            continue

        # assembly_build = col_list[3]  # MAF ASSEMBLY

        # chromosome is expected to be one or two number or one letter
        if hasattr(data, "chromosome"):
            chrom = data.chromosome  # MAF CHROM
        else:
            raise ValueError("Chromosome is not defined in MAF file")

        if not keep_filtered and not passes_filter(getattr(data, "filter", None)):
            n_filtered += 1
            continue

        sample = resolve_sample(data)
        x, y = resolve_alleles(data)  # MAF REF and variant allele

        if y is None:
            continue
        # skip if found unexpected nucleotide characters
        if len(set([x, y]) - set(nucleotides)) > 0:
            continue

        if hasattr(data, "start_position"):
            try:
                pos = int(data.start_position)
            except ValueError:
                raise ValueError("Start position is not a number in MAF file")
        else:
            raise ValueError("Start position is not defined in MAF file")

        # assuming that reference strand for reported mutations is always '+'!

        # transcript strand could be anything
        if hasattr(data, "transcript_strand"):
            try:
                transcript_strand = _normalize_transcript_strand(data.transcript_strand)
            except ValueError as e:
                # unusable strand in this row: skip it and keep reading the file
                skipped_rows.append((tracker.line_number, str(e)))
                N_skipped += 1
                continue
        else:
            # this is an incorrect assumption about transcript strand
            transcript_strand = "+"

        raw_mutations[sample].append((chrom, pos, transcript_strand, x, y))

    report_malformed_rows(skipped_rows, "MAF")
    report_filtered(n_filtered, "MAF")

    mutations, mutations_with_context, n_skipped, context_stats = _assemble_mutations(
        raw_mutations, asm, window_size
    )
    N_skipped += n_skipped

    N_loaded = 0
    for sample, sample_mutations in mutations.items():
        N_loaded += int(sum(sample_mutations.values()))

    processing_stats = {
        "loaded": N_loaded,
        "skipped": N_skipped,
        "nsamples": len(mutations.keys()),
        "format": "MAF",
        **context_stats,
    }

    return mutations, mutations_with_context, processing_stats


def read_VCF_with_context_window(infile, asm, window_size, keep_filtered=False):
    """
    Read VCF file and extract context of mutations for assembly asm and window +/- window_size around each mutation
    returns mutations, mutations_with_context, processing_stats
    """
    raw_mutations = defaultdict(list)

    N_skipped = 0
    n_filtered = 0
    # N_skipped_indels = 0

    sample = "VCF"

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

        # FILTER is the seventh column. A row that stops before it has
        # recorded no filters, so there is nothing to honour.
        if not keep_filtered and len(col_list) > 6 and not passes_filter(col_list[6]):
            n_filtered += 1
            continue

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

        transcript_strand = "+"
        raw_mutations[sample].append((chrom, pos, transcript_strand, x, y))
    # print("RAW", raw_mutations)
    # print("INDELS", N_skipped)

    report_filtered(n_filtered, "VCF")

    mutations, mutations_with_context, n_skipped, context_stats = _assemble_mutations(
        raw_mutations, asm, window_size
    )
    N_skipped += n_skipped

    N_loaded = 0
    for sample, sample_mutations in mutations.items():
        N_loaded += int(sum(sample_mutations.values()))

    nsamples = len(mutations.keys())
    processing_stats = {
        "loaded": N_loaded,
        "skipped": N_skipped,
        "format": "VCF",
        "nsamples": nsamples,
        **context_stats,
    }
    return mutations, mutations_with_context, processing_stats
