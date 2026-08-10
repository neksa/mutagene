import logging
import tarfile
from collections import defaultdict

from mutagene.io.profile import read_profile_str

logger = logging.getLogger(__name__)


def read_cohort_size_from_profile_file(profile_file):
    try:
        with open(profile_file) as f:
            profile_str = f.read()
            return read_cohort_size_from_profile_str(profile_str)
    except OSError:
        logger.warning("Could not read profile file")
        return 0


def read_cohort_size_from_profile_str(profile_str):
    for line in profile_str.splitlines():
        if not line.startswith("#"):
            continue
        fields = line.strip().split()
        # a comment line can hold anything, only #NSAMPLES <n> is meaningful here
        if len(fields) != 2 or fields[0] != "#NSAMPLES":
            continue
        try:
            return int(fields[1])
        except ValueError:
            # some published cohorts carry '#NSAMPLES None' and an unknown
            # cohort size must not take down the whole run
            logger.warning(f"Cohort size is not a number, treating it as unknown: {fields[1]}")
            return 0
    return 0


def list_cohorts_in_tar(tar_fname):
    """Returns a multiline string formatted list of cohorts contained in tar file"""
    suffix = ".aa_mutations.txt"
    cohorts = set()
    with tarfile.open(tar_fname, "r:*") as tar:
        for t in tar:
            name = t.name.rsplit("/", 1)[-1]
            # skip AppleDouble sidecars and other hidden entries
            if name.startswith("."):
                continue
            if name.lower().endswith(suffix):
                cohorts.add(name[: -len(suffix)])
    return "\n".join("\t" + cohort for cohort in sorted(cohorts))


def read_aa_mutations_map(aa_str):
    mutations = defaultdict(dict)
    for line in aa_str.splitlines():
        if len(line) == 0:
            continue
        fields = line.split()
        if len(fields) != 3:
            continue
        gene, mut, count = fields
        count = int(count)
        mutations[gene][mut] = count
    return mutations


def read_na_mutations_map(na_str):
    mutations = defaultdict(dict)
    for line in na_str.splitlines():
        if len(line) == 0:
            continue
        fields = line.split()
        if len(fields) != 5:
            continue
        chrom, pos, ref, alt, count = fields
        count = int(count)
        mutations[chrom][(pos, ref, alt)] = count
    return mutations


def read_cohort_mutations_from_tar(tar_fname, cohort):
    """Loads up profile, cohort size, aa mutations and na mutations from precalculated cohorts tar file"""
    aa_mutations = {}
    na_mutations = {}
    profile = []
    cohort_size = 0
    prefix = f"{cohort.lower()}."
    with tarfile.open(tar_fname, "r:*") as tar:
        for t in tar:
            name = t.name.rsplit("/", 1)[-1]
            # skip AppleDouble sidecars and other hidden entries
            if name.startswith("."):
                continue
            haystack = name.lower()
            # match on the file name so that bundles without a top level
            # directory load the same cohorts that list_cohorts_in_tar shows
            if haystack.startswith(prefix):
                if haystack.endswith(".profile"):
                    profile_str = tar.extractfile(t).read().decode("utf-8")
                    profile = read_profile_str(profile_str)
                    cohort_size = read_cohort_size_from_profile_str(profile_str)
                if haystack.endswith(".aa_mutations.txt"):
                    aa_mutations_str = tar.extractfile(t).read().decode("utf-8")
                    aa_mutations = read_aa_mutations_map(aa_mutations_str)
                if haystack.endswith(".dna_mutations.txt"):
                    na_mutations_str = tar.extractfile(t).read().decode("utf-8")
                    na_mutations = read_na_mutations_map(na_mutations_str)
    return profile, cohort_size, aa_mutations, na_mutations
