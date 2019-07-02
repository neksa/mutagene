import tarfile
from collections import defaultdict
from mutagene.io.profile import read_profile_str

import logging
logger = logging.getLogger(__name__)


def read_cohort_size_from_profile_file(profile_file):
    try:
        with open(profile_file) as f:
            profile_str = f.read()
            return read_cohort_size_from_profile_str(profile_str)
    except IOError:
        logger.warning("Could not read profile file")
        return 0


def read_cohort_size_from_profile_str(profile_str):
    for line in profile_str.splitlines():
        if line.startswith("#"):
            # print(line)
            a, b = line.strip().split()
            if a == "#NSAMPLES":
                return int(b)
    return 0


def list_cohorts_in_tar(tar_fname):
    """ Returns a multiline string formatted list of cohorts contained in tar file """
    cohorts = []
    with tarfile.open(tar_fname, 'r:gz') as tar:
        for t in tar:
            haystack = t.name.lower()
            if haystack.endswith(".aa_mutations.txt"):
                cohorts.append("\t" + haystack.split("/")[1].split('.')[0])
    return "\n".join(cohorts)


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
    """ Loads up profile, cohort size, aa mutations and na mutations from precalculated cohorts tar file """
    aa_mutations = {}
    na_mutations = {}
    profile = []
    cohort_size = 0
    with tarfile.open(tar_fname, 'r:gz') as tar:
        for t in tar:
            haystack = t.name.lower()
            needle = "/{}.".format(cohort.lower())
            if haystack.find(needle) != -1:
                if haystack.endswith('.profile'):
                    profile_str = tar.extractfile(t).read().decode('utf-8')
                    profile = read_profile_str(profile_str)
                    cohort_size = read_cohort_size_from_profile_str(profile_str)
                if haystack.endswith('.aa_mutations.txt'):
                    aa_mutations_str = tar.extractfile(t).read().decode('utf-8')
                    aa_mutations = read_aa_mutations_map(aa_mutations_str)
                if haystack.endswith('.dna_mutations.txt'):
                    na_mutations_str = tar.extractfile(t).read().decode('utf-8')
                    na_mutations = read_na_mutations_map(na_mutations_str)
    return profile, cohort_size, aa_mutations, na_mutations
