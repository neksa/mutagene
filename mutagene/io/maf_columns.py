"""Column resolution shared by the MAF readers.

MutaGene has two MAF parsers -- ``read_MAF_with_context_window`` (used by rank,
motif and the webapp) and ``read_MAF_profile`` (used by profile). They grew
their own header and allele handling and drifted until they accepted different
files: issues #91, #92 and #98 were each reported as one command working on a
MAF that another command rejected. Resolving columns here keeps them in step, so
a file either parses for every subcommand or for none of them.
"""


def normalize_header(header):
    """Map MAF column names onto namedtuple-safe field names.

    Column names vary in case between MAF producers, and some carry dots
    (``Func.refGene``), which are not valid Python identifiers.
    """
    return tuple(name.lower().replace(".", "_") for name in header)


def resolve_sample(data):
    """Return the sample identifier for a MAF row."""
    for field in ("tumor_sample_barcode", "sample_id"):
        value = getattr(data, field, None)
        if value is not None:
            return value
    raise ValueError("Sample ID is not defined in MAF file")


def resolve_alleles(data):
    """Return ``(reference_allele, variant_allele)`` for a MAF row.

    ``Tumor_Seq_Allele1`` commonly repeats the reference allele, so the variant
    is whichever tumour allele differs from the reference, preferring allele 2.
    Either tumour column on its own is enough: requiring *both* is what made
    #98 reject files that the profile reader read without complaint.

    The variant is ``None`` when no tumour allele differs from the reference,
    which callers treat as a row carrying no substitution.
    """
    reference = getattr(data, "reference_allele", None)
    if reference is None:
        raise ValueError("Reference allele is not defined in MAF file")

    variant = getattr(data, "variant_allele", None)
    if variant is not None:
        return reference, variant

    allele1 = getattr(data, "tumor_seq_allele1", None)
    allele2 = getattr(data, "tumor_seq_allele2", None)
    if allele1 is None and allele2 is None:
        raise ValueError("Variant allele is not defined in MAF file")

    variant = allele1 if allele1 != reference else None
    if allele2 is not None and allele2 != reference:
        variant = allele2
    return reference, variant
