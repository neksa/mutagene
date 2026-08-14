"""Counters describing how well mutations matched the reference assembly.

Both context extractors report these, so a caller can tell "the wrong genome was
chosen" apart from ordinary data loss without scraping log messages. Counting
warning records was what made the webapp's mismatch detector silently useless
(#99): it watched a logger during a window in which the message it wanted was
never emitted, and the message it counted is not produced by the profile path at
all.

``ref_mismatches`` is the assembly signal. It counts only mutations whose
reference allele matches *neither* the assembly base nor its complement, which
is what a wrong assembly produces. A MAF that simply reports alleles on the
reverse strand is counted under ``reverse_strand_ref`` instead: it is normal,
handled, and must not be mistaken for a wrong genome.
"""

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# A wrong assembly puts the reference allele at odds with the sequence for a
# large share of mutations, not for a handful. Requiring a rate as well as a
# count keeps a few dozen genuinely odd rows in a multi-million-mutation file
# from being announced as the wrong genome.
MISMATCH_MIN_COUNT = 20
MISMATCH_MIN_RATE = 0.05

_ALTERNATIVE_ASSEMBLY = {"hg19": "hg38", "hg38": "hg19", "mm9": "mm10", "mm10": "mm9"}

CONTEXT_STAT_KEYS = (
    "ref_mismatches",
    "reverse_strand_ref",
    "chromosome_not_found",
    "read_errors",
)


def new_context_stats():
    """Return a zeroed set of context counters."""
    return dict.fromkeys(CONTEXT_STAT_KEYS, 0)


def merge_context_stats(total, part):
    """Add the counters in *part* into *total* in place and return it."""
    for key in CONTEXT_STAT_KEYS:
        total[key] = total.get(key, 0) + part.get(key, 0)
    return total


def mismatch_rate(stats, loaded):
    """Fraction of examined mutations whose REF matched neither strand.

    The denominator is every mutation that was looked at, mismatches included,
    so a file where nothing survived still reports a rate near 1 rather than
    dividing by zero.
    """
    examined = loaded + stats.get("ref_mismatches", 0)
    if examined <= 0:
        return 0.0
    return stats.get("ref_mismatches", 0) / examined


def assembly_name(genome):
    """Reduce a genome argument, which may be a full 2bit path, to its name."""
    name = Path(str(genome)).name
    return name[: -len(".2bit")] if name.endswith(".2bit") else name


def assembly_mismatch_warning(stats, loaded, genome):
    """Describe a suspected wrong genome assembly, or return None.

    The payload keeps the ``mismatch_count`` and ``message`` keys the webapp
    already stores and renders; the rest is additive.
    """
    count = stats.get("ref_mismatches", 0)
    rate = mismatch_rate(stats, loaded)
    if count <= MISMATCH_MIN_COUNT or rate < MISMATCH_MIN_RATE:
        return None

    name = assembly_name(genome)
    alternative = _ALTERNATIVE_ASSEMBLY.get(name)
    suggestion = f" Please try {alternative}." if alternative else ""
    return {
        "mismatch_count": count,
        "examined": loaded + count,
        "mismatch_rate": rate,
        "genome": name,
        "message": (
            f"Found {count} reference allele mismatches "
            f"({rate:.0%} of mutations examined). This suggests the wrong genome "
            f"assembly was selected.{suggestion}"
        ),
    }


def report_assembly_mismatch(stats, loaded, genome):
    """Log what the mismatch counters say and return the warning payload, if any."""
    count = stats.get("ref_mismatches", 0)
    if count == 0:
        return None

    logger.warning(f"{count} mutations have a reference allele that matches neither strand")
    warning = assembly_mismatch_warning(stats, loaded, genome)
    if warning is not None:
        logger.error(warning["message"])
    return warning
