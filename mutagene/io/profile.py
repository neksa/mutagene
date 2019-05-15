# from tqdm import tqdm
from collections import defaultdict
from mutagene.dna import nucleotides
from mutagene.dna import complementary_trinucleotide, complementary_nucleotide

import os
import numpy as np

import logging
logger = logging.getLogger(__name__)


def read_profile_file(profile_file):
    """Read profile from file

    Format: T[A>C]G frequency

    Arguments:
        profile_file {str} -- profile file name

    Returns:
        mutations, stats --
    """
    try:
        with open(profile_file) as f:
            profile_str = f.read()
            return read_profile_str(profile_str)
    except IOError:
        logger.warning("Could not read profile file " + profile_file)
        return None


def read_profile_str(profile_str):
    mutations = defaultdict(float)
    for line in profile_str.splitlines():
        if len(line) == 0:
            continue
        if line.startswith('#'):
            continue
        fields = line.strip().upper().split()
        if len(fields) != 2:
            return None, None

        if len(fields[0]) != 7:
            return None, None

        if fields[0][1] != "[" or fields[0][3] != ">" or fields[0][5] != "]":
            return None, None

        p5, _, x, _, y, _, p3 = tuple(fields[0])

        if p5 not in nucleotides or p3 not in nucleotides or y not in nucleotides:
            return None, None

        if x not in "TC":
            return None, None

        try:
            f = float(fields[1])
        except:
            return None, None

        mutations[p5 + p3 + x + y] = f

    values = []
    for p5 in nucleotides:
        for p3 in nucleotides:
            for x in "CT":
                for y in nucleotides:
                    if x != y:
                        values.append(mutations.get(p5 + p3 + x + y, 0.0))
    return values


def read_signatures(n_signatures):
    signatures_dict = {5: 'A', 10: 'B', 30: 'C'}
    assert n_signatures in signatures_dict

    dirname = os.path.dirname(os.path.realpath(__file__))

    W = []
    signature_names = []
    for i in range(n_signatures):
        fname = dirname + "/../data/signatures/{}_{}.profile".format(signatures_dict[n_signatures], i + 1)
        profile = read_profile_file(fname)
        W.append(profile)
        signature_names.append("{}".format(i + 1))

    W = np.array(W).T
    return W, signature_names


def write_profile(profile_file, p, counts=True):
    with open(profile_file, 'w') as o:
        write_profile_file(o, p, counts)


# def format_profile_dict(values):
#     attrib = get_signature_attributes_dict()
#     result = {}
#     for i, v in enumerate(values):
#         x, y = attrib[i]['mutation']
#         p5, p3 = attrib[i]['context']
#         result["{}[{}>{}]{}".format(p5, x, y, p3)] = v
#     return result


def format_profile(values, counts=False):
    if counts:
        format_str = "{}[{}>{}]{}\t{:d}\n"
    else:
        format_str = "{}[{}>{}]{}\t{:.12f}\n"
    attrib = get_profile_attributes_dict()
    result = ""
    for i, v in enumerate(values):
        x, y = attrib[i]['mutation']
        p5, p3 = attrib[i]['context']
        if counts:
            v = int(v)
        result += format_str.format(p5, x, y, p3, v)
    return result


def write_profile_file(file_handle, p, counts=True):
    formatted_profile = format_profile(p, counts)
    file_handle.write(formatted_profile)


def get_profile_attributes_dict(signature_order=False):
    attribs = []
    if signature_order:
        for x in "CT":
            for y in nucleotides:
                if x != y:
                    for p5 in nucleotides:
                        for p3 in nucleotides:
                            attribs.append({
                                'mutation': x + y,
                                'context': p5 + p3
                            })
    else:
        for p5 in nucleotides:
            for p3 in nucleotides:
                for x in "CT":
                    for y in nucleotides:
                        if x != y:
                            attribs.append({
                                'mutation': x + y,
                                'context': p5 + p3
                            })
    return attribs


def get_attributes():
    attribs = get_profile_attributes_dict()
    for i, v in enumerate(attribs):
        x, y = v['mutation']
        p5, p3 = v['context']
        attribs[i]['mut'] = "{}{}".format(x, y),
        attribs[i]['mutation'] = "{}{}{} → {}{}{}".format(p5, x, p3, p5, y, p3),
        attribs[i]['context'] = "{}{}".format(p5, p3)
    return attribs


def profile_to_dict(counts_profile):
    mutation_prob = defaultdict(dict)
    total_count = sum(counts_profile)
    freqs = [x / total_count for x in counts_profile]

    for i, attrib in enumerate(get_profile_attributes_dict()):
        x, y = attrib['mutation']
        p5, p3 = attrib['context']
        tri = p5 + x + p3
        mutation_prob[tri][y] = freqs[i]
        mutation_prob[complementary_trinucleotide[tri]][complementary_nucleotide[y]] = freqs[i]
    return mutation_prob
