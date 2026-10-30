"""Validate the VCF sample layout before reading variants."""


def validate_vcf_header(line):
    """Reject multiple sample columns until the parser supports sample selection."""
    if line.startswith("#CHROM") and len(line.rstrip().split("\t")) > 10:
        raise ValueError(
            "Multiple VCF samples are not supported. "
            "Select one tumor sample and retain its variant calls before analysis."
        )


def has_reference_genotype(fields):
    """Return whether the single sample explicitly contains only reference alleles."""
    if len(fields) != 10:
        return False
    keys = fields[8].split(":")
    if "GT" not in keys:
        return False
    values = fields[9].split(":")
    index = keys.index("GT")
    if index >= len(values):
        return False
    alleles = values[index].replace("|", "/").split("/")
    return all(allele == "0" for allele in alleles)
