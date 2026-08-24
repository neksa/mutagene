"""Honouring the FILTER field that variant callers write.

FILTER is a VCF column: the caller writes ``PASS`` for a variant it stands
behind, ``.`` when it applied no filters, and otherwise the names of the filters
the variant failed (``germline``, ``weak_evidence``, ``clustered_events`` and so
on, semicolon-separated). GDC-style MAFs carry the same field forward in a
column named FILTER, so this is not a MAF-specific idea and both readers use it.

Reading a caller's rejected variants as though they were accepted silently
changes what an analysis is about. A Mutect2 file in which 4,637 of 4,675 rows
are marked ``germline`` produces a "mutational profile" that is mostly germline
variation, and the decomposition duly reports germline-contamination signatures.

Only a column named exactly FILTER is used. GDC's separate ``GDC_FILTER`` is
deliberately ignored: it mixes quality flags with annotations such as
``NonExonic``, and non-exonic mutations are perfectly good input to a mutational
profile. Dropping them would lose real data to a filter that never meant that.
"""

import logging

logger = logging.getLogger(__name__)

# "PASS" is an accepted variant; "." and empty mean no filters were applied.
_PASSING_VALUES = frozenset({"PASS", ".", ""})


def add_filter_argument(parser):
    """Attach --keep-filtered to a subcommand that reads variants."""
    parser.add_argument(
        "--keep-filtered",
        action="store_true",
        help="Include variants the caller rejected in its FILTER column "
        "(germline, weak_evidence and so on). By default only PASS variants "
        "and files with no FILTER column are used",
    )


def passes_filter(value):
    """True when a FILTER value describes a variant the caller accepted.

    A missing column gives None, which passes: a file that records no filters
    has not rejected anything, and refusing all of it would be wrong.
    """
    if value is None:
        return True
    return value.strip() in _PASSING_VALUES


def report_filtered(n_filtered, file_format):
    """Say how many variants the caller's own FILTER field excluded."""
    if not n_filtered:
        return
    logger.info(
        f"Excluded {n_filtered} variants marked as filtered in the {file_format} "
        "file's FILTER column. Use --keep-filtered to include them"
    )
