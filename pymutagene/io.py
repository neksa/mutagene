import numpy as np
from collections import defaultdict

from .dna import nucleotides


def read_signatures(n_signatures):
    signatures_dict = {5: 'A', 10: 'B', 30: 'C'}
    assert n_signatures in signatures_dict

    W = []
    signature_names = []
    for i in range(n_signatures):
        fname = "data/signatures/{}_{}.profile".format(signatures_dict[n_signatures], i + 1)
        profile = read_profile(fname)
        W.append(profile)
        signature_names.append("{}".format(i + 1))

    W = np.array(W).T
    return W, signature_names


def write_profile(profile_file, p, counts=True):
    with open(profile_file, 'w') as o:
        formatted_profile = format_profile(p, counts)
        o.write(formatted_profile)


def read_profile(profile_file):
    """Read profile from file

    Format: T[A>C]G frequency

    Arguments:
        profile_file {str} -- profile file name

    Returns:
        mutations, stats --
    """
    try:
        with open(profile_file) as f:
            mutations = defaultdict(float)
            for line in f:
                if len(line) == 0:
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
                # print(p5, p3, x, y)

            values = []
            for p5 in nucleotides:
                for p3 in nucleotides:
                    for x in "CT":
                        for y in nucleotides:
                            if x != y:
                                values.append(mutations.get(p5 + p3 + x + y, 0.0))
            return values
    except IOError:
        return None


def get_dummy_signatures_lists():
    """
    Generate 6 dummy signatures
    Each will have uniform non-zero frequencies corresponding to one mutation type
    Format them as lists
    """
    dummy_signatures = []
    for mutation in (("C", "A"), ("C", "T"), ("C", "G"), ("T", "A"), ("T", "C"), ("T", "G")):
        values = []
        for p5 in nucleotides:
            for p3 in nucleotides:
                for x in "CT":
                    for y in nucleotides:
                        if x != y:
                            if mutation == (x, y):
                                values.append(1.0 / 16.0)
                            else:
                                values.append(0.0)
        name = mutation[0] + " to " + mutation[1]
        dummy_signatures.append((name, values))
    return dummy_signatures


def get_signature_attributes_dict(signature_order=False):
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
    attribs = get_signature_attributes_dict()
    for i, v in enumerate(attribs):
        x, y = v['mutation']
        p5, p3 = v['context']
        attribs[i]['mut'] = "{}{}".format(x, y),
        attribs[i]['mutation'] = "{}{}{} → {}{}{}".format(p5, x, p3, p5, y, p3),
        attribs[i]['context'] = "{}{}".format(p5, p3)
    return attribs


def format_profile(values, counts=False):
    if counts:
        format_str = "{}[{}>{}]{}\t{:d}\n"
    else:
        format_str = "{}[{}>{}]{}\t{:.12f}\n"
    attrib = get_signature_attributes_dict()
    result = ""
    for i, v in enumerate(values):
        x, y = attrib[i]['mutation']
        p5, p3 = attrib[i]['context']
        if counts:
            v = int(v)
        result += format_str.format(p5, x, y, p3, v)
    return result


# def format_profile_dict(values):
#     attrib = get_signature_attributes_dict()
#     result = {}
#     for i, v in enumerate(values):
#         x, y = attrib[i]['mutation']
#         p5, p3 = attrib[i]['context']
#         result["{}[{}>{}]{}".format(p5, x, y, p3)] = v
#     return result
