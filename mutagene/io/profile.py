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
def read_COSMICv3_signatures():
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


def read_KUCAB_signatures():
    """
    A compendium of mutational signatures of environmental agents Kucab et al.  Serena Nik-Zainal
    DOI: 10.17632/m7r4msjb4c.2
    https://data.mendeley.com/datasets/m7r4msjb4c/2#folder-9f2942fa-935b-4ac3-aa26-35784fb68a85

    MutationType    Potassium bromate (875 uM)  DBADE (0.109 uM)    Formaldehyde (120 uM)   Semustine (150 uM)  Temozolomide (200 uM)   DMH (11.6 mM) + S9  Benzidine (200 uM)  DBP (0.0039 uM) MX (7 uM) + S9  Methyleugenol (1.25 mM) 4-ABP (300 uM) + S9 DBPDE (0.000156 uM) DBP (0.0313 uM) + S9    DBADE (0.0313 uM)   1,8-DNP (0.125 uM)  BPDE (0.125 uM) MNU (350 uM)    ENU (400 uM)    Cyclophosphamide (18.75 uM) + S9    BaP (0.39 uM) + S9  6-Nitrochrysene (12.5 uM) + S9  AAI (1.25 uM)   Potassium bromate (260 uM)  6-Nitrochrysene (0.78 uM)   Ellipticine (0.375 uM) + S9 DBA (75 uM) + S9    PhIP (3 uM) + S9    AFB1 (0.25 uM) + S9 3-NBA (0.025 uM)    1,6-DNP (0.09 uM)   5-Methylchrysene (1.6 uM) + S9  Furan (100 mM) + S9 SSR (1.25 J)    AAII (37.5 uM)  Propylene oxide (10 mM) N-Nitrosopyrrolidine (50 mM)    Mechlorethamine (0.3 uM)    DES (0.938 mM)  DMS (0.078 mM)  Cisplatin (3.125 uM)    OTA (0.08 uM) + S9  Carboplatin (5 uM)  DBAC (5 uM) + S9    Temozolomide (200 uM).1 Cisplatin (12.5 uM) AZD7762 (1.625 uM)  3-NBA (0.1 uM)  PhIP (4 uM) + S9    BaP (2 uM) + S9 6-Nitrochrysene (50 uM) + S9    6-Nitrochrysene (50 uM) 1,8-DNP (8 uM)  DBPDE (0.000625 uM) Control
    A[C>A]A 0.109398510425891   0.0240330588511261  0.0181865586268237  0.000122421190463747    0   0.00665025743214769 0.0564462237667514  0   0.0494977586860332  0.111458492403729   0.00483165426221388 0.0002218019424568  0.00869925685873544 0.0264652966450285  0.0195599544767651  0.0237017339519746  0.00262256704780487 0.0130480902574389  0.00555898815479868 0.0406568174086766  0.0100871864208784  0.00428520931740606 0.114738757387778   0.00515014598157648 0.0482801050571615  0.0323606741641645  0.0279971115871403  0.0379396506608179  0.000340069771278152    0.000751844456779883    0.0286278352872072  0.0503953181060661  2.22941465639634e-05    0.000162365803469566    0.0104143358776828  0.00685630323850839 0.000133921576621699    0.014166642071802   0.000461507809128427    0.00279492689318811 0.11584205195812    0.00194170922837543 0.00673772703708239 0.000256609057933161    1.12145020078519e-05    0.000395123818576932    0.0258298926507019  0.032667880519281   0.0283815533748292  0.00855231498108648 0.00555769823738427 0.0180349250334742  0.0186419336589614  0.0605073306958345
    """
    mutations = defaultdict(dict)
    dirname = os.path.dirname(os.path.realpath(__file__))
    fname = os.path.normpath(dirname + "/../data/signatures/Kucab.txt")
    signature_names = []
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            fields = line.strip().split("\t")
            if i == 0:
                signature_names = fields[1:]
                continue

            # parse A[C>A]A
            m = fields[0]
            p5 = m[0]
            assert m[1] == '['
            x = m[2]
            assert m[3] == '>'
            y = m[4]
            assert m[5] == ']'
            p3 = m[6]

            for j, val in enumerate(fields[1:]):
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


def _read_mutagene_signatures(name, n_signatures):
    """ Read signatures preprocessed for mutagene website / package """
    dirname = os.path.dirname(os.path.realpath(__file__))
    W = []
    signature_names = []
    for i in range(n_signatures):
        fname = dirname + "/../data/signatures/{}_{}.profile".format(name, i + 1)
        profile = read_profile_file(fname)
        W.append(profile)
        signature_names.append("{}".format(i + 1))

    W = np.array(W).T
    return W, signature_names


def read_MGA_signatures():
    return _read_mutagene_signatures('A', 5)


def read_MGB_signatures():
    return _read_mutagene_signatures('B', 10)


def read_COSMICv2_signatures():
    return _read_mutagene_signatures('C', 30)


def read_signatures(name, only=None):
    """
    Retrieve a set of signatures by its name or number of signatures specified as 'str'.
    Returns tuple: numpy matrix and list of names
    """

    # number of sigatures matching to names
    signatures_dict = {
        '5': 'MGA',
        '10': 'MGB',
        '30': 'COSMICv2',
        '49': 'COSMICv3',
        '53': 'KUKAB'
    }
    inv_signatures_dict = dict(zip(signatures_dict.values(), signatures_dict.keys()))

    if name in signatures_dict:
        name = signatures_dict[name]

    assert name in inv_signatures_dict, "Unknown name for a signature set: {}".format(name)

    function_name = "read_{}_signatures".format(name)
    W, signature_names = globals()[function_name]()

    # filter by list of 'only' signatures,
    # append signature set prefix to signature names
    delete_signatures = []
    filtered_signature_names = []
    for i, x in enumerate(signature_names):
        if only is not None:
            if x not in only:
                delete_signatures.append(i)
                continue
        filtered_signature_names.append("{}-{}".format(name, x))

    if only is not None:
        W = np.delete(W, np.array(delete_signatures), axis=1)

    return W, filtered_signature_names


def write_profile(profile_file, p, counts=True):
    with open(profile_file, 'w') as o:
        write_profile_file(o, p, counts)


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
