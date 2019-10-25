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


"""
Type,SubType,SBS1,SBS2,SBS3,SBS4,SBS5,SBS6,SBS7a,SBS7b,SBS7c,SBS7d,SBS8,SBS9,SBS10a,SBS10b,SBS11,SBS12,SBS13,SBS14,SBS15,SBS16,SBS17a,SBS17b,SBS18,SBS19,SBS20,SBS21,SBS22,SBS23,SBS24,SBS25,SBS26,SBS27,SBS28,SBS29,SBS30,SBS31,SBS32,SBS33,SBS34,SBS35,SBS36,SBS37,SBS38,SBS39,SBS40,SBS41,SBS42,SBS43,SBS44,SBS45,SBS46,SBS47,SBS48,SBS49,SBS50,SBS51,SBS52,SBS53,SBS54,SBS55,SBS56,SBS57,SBS58,SBS59,SBS60,SBS84,SBS85
C>A,ACA,8.86E-04,5.80E-07,2.08E-02,4.22E-02,1.20E-02,4.25E-04,6.70E-05,2.33E-03,4.83E-03,4.04E-05,4.41E-02,5.58E-04,2.19E-03,1.82E-04,1.46E-04,4.52E-03,1.82E-03,1.12E-03,9.44E-04,1.60E-02,2.07E-03,6.08E-04,5.15E-02,1.27E-03,6.19E-04,1.57E-04,6.01E-03,8.35E-04,3.64E-02,9.90E-03,8.73E-04,5.21E-03,7.84E-04,6.32E-02,1.80E-03,9.54E-03,2.23E-02,3.11E-03,4.87E-03,8.83E-03,2.52E-02,3.95E-03,1.28E-02,1.17E-02,2.82E-02,2.11E-03,1.16E-03,2.95E-02,7.68E-18,9.11E-03,4.40E-03,6.78E-02,8.55E-04,2.51E-02,1.19E-01,1.41E-01,1.52E-02,5.38E-03,2.16E-03,5.88E-03,1.26E-02,1.23E-02,5.89E-02,3.59E-03,6.15E-03,0.003471994832,0.006080257390
"""
"""
Sequencing artifacts:

Possible sequencing artefacts
"""
def read_COSMIC3_signatures():
    sequencing_artifacts = [
        "SBS27", "SBS43", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49", "SBS50", "SBS51",
        "SBS52", "SBS53", "SBS54", "SBS55", "SBS56", "SBS57", "SBS58", "SBS59", "SBS60"]

    mutations = defaultdict(dict)
    dirname = os.path.dirname(os.path.realpath(__file__))
    fname = dirname + "/../data/signatures/sigProfiler_SBS_signatures_2019_05_22.csv"
    signature_names = []
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            fields = line.strip().split(',')   
            if i == 0:
                signature_names = fields[2:]
                continue

            x, y = fields[0].split(">")
            p5, _, p3 = tuple(fields[1])

            for j, val in enumerate(fields[2:]):
                val = float(val)
                mutations[j][p5 + p3 + x + y] = val

    W = []
    for j in range(len(signature_names)):
        profile = []
        for p5 in nucleotides:
            for p3 in nucleotides:
                for x in "CT":
                    for y in nucleotides:
                        if x != y:
                            profile.append(mutations[j].get(p5 + p3 + x + y, 0.0))
        W.append(profile)
    W = np.array(W).T
    return W, signature_names


def read_signatures(n_signatures):
    signatures_dict = {5: 'A', 10: 'B', 30: 'C', 49: 'C3'}
    assert n_signatures in signatures_dict

    if n_signatures == 49:
        return read_COSMIC3_signatures()

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
