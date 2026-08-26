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

By default only a column named exactly FILTER is used. GDC's separate
``GDC_FILTER`` is not read: it mixes quality flags with annotations such as
``NonExonic``, and non-exonic mutations are perfectly good input to a mutational
profile. Dropping them would lose real data to a filter that never meant that.

Some converters do not produce a FILTER column at all and leave the caller's
verdict in a passthrough column under a name that carries no meaning -- ANNOVAR
writes Mutect2's verdicts into ``Otherinfo10``, for instance. ``--filter-column``
names the column to read in that case. It is deliberately not guesswork: a
column requested and not found is an error, because silently reading nothing is
the failure this whole module exists to prevent.
"""

import logging

logger = logging.getLogger(__name__)

# "PASS" is an accepted variant; "." and empty mean no filters were applied.
_PASSING_VALUES = frozenset({"PASS", ".", ""})


DEFAULT_FILTER_COLUMN = "FILTER"


def add_filter_argument(parser):
    """Attach the FILTER options to a subcommand that reads variants."""
    parser.add_argument(
        "--filter-column",
        metavar="NAME",
        type=str,
        default=DEFAULT_FILTER_COLUMN,
        help="Name of the column holding the caller's FILTER verdict in a MAF "
        "file, for converters that write it elsewhere (ANNOVAR uses "
        "Otherinfo10). Default: FILTER",
    )
    parser.add_argument(
        "--keep-filtered",
        action="store_true",
        help="Include variants the caller rejected in its FILTER column "
        "(germline, weak_evidence and so on). By default only PASS variants "
        "and files with no FILTER column are used",
    )


def field_name(column):
    """The namedtuple field a MAF column name becomes after header normalisation."""
    return column.lower().replace(".", "_")


def check_filter_column(header, column, keep_filtered=False):
    """Fail early when a specifically requested FILTER column is not there.

    An absent default FILTER column is ordinary -- most MAFs have none, and the
    file has then rejected nothing. A column the caller *named* is different: if
    it is missing, nothing would be filtered, and a run that silently keeps every
    germline variant is exactly the outcome being guarded against.
    """
    if keep_filtered or column == DEFAULT_FILTER_COLUMN:
        return
    if field_name(column) not in header:
        raise ValueError(
            f"No column named {column!r} in this file, so --filter-column has "
            "nothing to read. Check the header, or drop the option to use FILTER"
        )


def filter_value(data, column):
    """The FILTER verdict for one row, or None when the column is absent."""
    return getattr(data, field_name(column), None)


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
