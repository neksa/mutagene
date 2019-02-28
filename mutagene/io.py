from collections import defaultdict, namedtuple
import json
import multiprocessing
import requests
import tarfile
import os
import csv
import urllib

import numpy as np
import twobitreader as tbr

from .dna import nucleotides, complementary_nucleotide

# TWOBIT_GENOMES_PATH = '/Users/agoncear/data/'


def read_signatures(n_signatures):
    signatures_dict = {5: 'A', 10: 'B', 30: 'C'}
    assert n_signatures in signatures_dict

    dirname = os.path.dirname(os.path.realpath(__file__))

    W = []
    signature_names = []
    for i in range(n_signatures):
        # fname = dirname + "/data/signatures/{}_{}.profile".format(signatures_dict[n_signatures], i + 1)
        fname = dirname + "/../data/signatures/{}_{}.profile".format(signatures_dict[n_signatures], i + 1)
        # print(fname)
        profile = read_profile_file(fname)
        W.append(profile)
        signature_names.append("{}".format(i + 1))

    W = np.array(W).T
    return W, signature_names


def write_profile(profile_file, p, counts=True):
    with open(profile_file, 'w') as o:
        write_profile_file(o, p, counts)


def write_profile_file(file_handle, p, counts=True):
    formatted_profile = format_profile(p, counts)
    file_handle.write(formatted_profile)


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
        # print(p5, p3, x, y)

    values = []
    for p5 in nucleotides:
        for p3 in nucleotides:
            for x in "CT":
                for y in nucleotides:
                    if x != y:
                        values.append(mutations.get(p5 + p3 + x + y, 0.0))
    return values


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
        return None


def read_cohort_size_from_profile_str(profile_str):
    for line in profile_str.splitlines():
        if line.startswith("#"):
            # print(line)
            a, b = line.strip().split()
            if a == "#NSAMPLES":
                return int(b)


def get_mutational_profile(mutations, counts=False):
    attrib = get_signature_attributes_dict()
    values = []
    total_mut_number = sum(mutations.values())
    for i, attr in enumerate(attrib):
        number = mutations[attr['context'] + attr['mutation']]
        # freq = 0.000001 * number / total_mut_number
        if counts:
            freq = number
        else:
            freq = number / float(total_mut_number)
        # trinucleotide = attr['context'][0] + attr['mutation'][0] + attr['context'][1]
        # trinucleotide_freq = exome_trinucleotide_freq[trinucleotide]
        # values.append(3.0 * freq / trinucleotide_freq)
        values.append(freq)
    # print(values)
    return values


def write_decomposition(fname, results, signature_ids):
    if type(results) is np.ndarray:
        h = results
    else:
        exposure_dict = {x['name']: x['score'] for x in results}
        exposure = [exposure_dict[name] for name in signature_ids]
        h = np.array(exposure)

    # FIXME: code duplication
    if isinstance(fname, (str, bytes)):
        with open(fname, 'w') as o:
            for i in range(h.shape[0]):
                o.write("{}\t{:.4f}\n".format(signature_ids[i], h[i]))
    else:
        for i in range(h.shape[0]):
            fname.write("{}\t{:.4f}\n".format(signature_ids[i], h[i]))


def read_decomposition(fname):
    signature_ids = []
    h = []

    try:
        with open(fname) as f:
            for line in f:
                a, b = line.strip().split()
                signature_ids.append(a)
                h.append(float(b))
    except:
        return None, None
    # except FileNotFoundError:
    #     return None, None

    return np.array(h), signature_ids


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


def get_context_twobit(mutations, twobit_file):
    """
    User twobitreader to get context of mutations
    """
    contexts = {}

    f = tbr.TwoBitFile(twobit_file)

    cn = complementary_nucleotide
    for (chrom, pos, x, y) in mutations:
        start = int(pos) - 1  # zero-based numbering
        chrom = str(chrom)
        chromosome = chrom if chrom.startswith('chr') else 'chr' + chrom

        nuc5 = 'N'
        nuc3 = 'N'
        nuc = 'N'
        if chromosome in f:
            try:
                seq = f[chromosome][start - 1: start + 2]  # +/- 1 nucleotide
                nuc5, nuc, nuc3 = tuple(seq.upper())
            except:
                nuc = 'N'
                # print(chromosome, x, nuc5, nuc, nuc3)
            if nuc != 'N' and nuc != x:
                if cn[nuc] == x:
                    nuc3 = cn[nuc5]
                    nuc5 = cn[nuc3]
                else:
                    nuc3 = nuc5 = 'N'
        else:
            # print("NO CHROM", chromosome)
            pass
        contexts[(chrom, pos)] = (nuc5, nuc3)
    return contexts


def get_context_batch(mutations, assembly, method='twobit'):
    """
        Get context for a list of mutations [(chrom, pos, x, y) ] format
    """
    if assembly is None:
        assembly = 38

    if method is None:
        method = 'twobit'

    methods = {
        'ensembl': get_context_ensembl,
        'twobit': get_context_twobit
    }

    contexts = methods[method](mutations, assembly)
    return contexts


def read_mutations(muts, fmt=None, asm=None):
    mutations = None
    processing_stats = None
    if fmt is not None:
        fmt = fmt.upper()

    if fmt is None or fmt == "AUTO" or fmt == "":
        # for line in muts.split("\n"):
        for line in muts:
            if line.startswith("#version 2."):
                fmt = "MAF"
                break
            if len(line.strip()) == 0 or line.startswith("#"):
                continue
            tabs = line.split()
            if len(tabs) == 2 and "[" in tabs[0] and "]" in tabs[0]:
                fmt = "PROFILE"
                break
            if (len(tabs) == 3 and tabs[1] == ">") or (len(tabs) == 1 and ">" in tabs[0]):
                fmt = "TRI"
                break
            if len(tabs) > 3:
                fmt = "VCF"
                if tabs[0].lower().startswith("chr"):
                    break  # yes, it's VCF
                if tabs[4].lower().startswith("chr"):
                    fmt = "MAF"
                    break
            continue

    print("DATA FORMAT:", fmt)

    if fmt is None or fmt == "AUTO" or fmt == "" or fmt == "TRI":
        try:
            mutations, processing_stats = read_trinucleotide(muts)
        except:
            pass

    if fmt == "VCF":
        mutations, processing_stats = read_VCF_profile(muts, asm)
    if fmt == "MAF":
        mutations, processing_stats = read_MAF_profile(muts, asm)
    # if fmt == "PROFILE":
    #     mutations, processing_stats = read_PROFILE(muts)

    return mutations, processing_stats


def read_MAF_profile(muts, asm=None):
    cn = complementary_nucleotide
    mutations = defaultdict(float)
    N_skipped = 0

    raw_mutations = []
    for i, line in enumerate(muts):
        if line.startswith("#"):
            continue
        if line.startswith("Hugo_Symbol"):
            continue
        if len(line) < 10:
            continue

        col_list = line.split()
        if len(col_list) < 13:
            continue

        try:
            # assembly_build = col_list[3]  # MAF ASSEMBLY
            # strand = col_list[7]    # MAF STRAND
            # ID = col_list[2]

            # chromosome is expected to be one or two number or one letter
            chrom = col_list[4]  # MAF CHROM
            if chrom.lower().startswith("chr"):
                chrom = chrom[3:]
            # if len(chrom) == 2 and chrom[1] not in "0123456789":
            #     chrom = chrom[0]

            pos = int(col_list[5])      # MAF POS START
            pos_end = int(col_list[6])  # MAF POS END
            x = col_list[10]            # MAF REF
            y1 = col_list[11]           # MAF ALT1
            y2 = col_list[12]           # MAF ALT2
        except:
            # raise
            continue

        if pos != pos_end:
            continue

        # skip if found unexpected nucleotide characters
        if len(set([x, y1, y2]) - set(nucleotides)) > 0:
            continue

        y = y1 if y1 != x else None
        y = y2 if y2 != x else y
        if y is None:
            continue

        raw_mutations.append((chrom, pos, x, y))

    MAX = 100000000
    if len(raw_mutations) > 0:
        if len(raw_mutations) > MAX:
            raw_mutations = raw_mutations[:MAX]

        contexts = get_context_batch(raw_mutations, asm)

        if contexts is None:
            return None, None

        if len(contexts) == 0:
            return None, None

        for (chrom, pos, x, y) in raw_mutations:
            p5, p3 = contexts.get((chrom, pos), ("N", "N"))

            if len(set([p5, x, y, p3]) - set(nucleotides)) > 0:
                # print("Skipping invalid nucleotides")
                N_skipped += 1
                continue

            if x in "CT":
                mutations[p5 + p3 + x + y] += 1.0
            else:
                # complementary mutation
                mutations[cn[p3] + cn[p5] + cn[x] + cn[y]] += 1.0

    N_loaded = int(sum(mutations.values()))
    processing_stats = {'loaded': N_loaded, 'skipped': N_skipped, 'format': 'MAF'}
    return mutations, processing_stats


def read_VCF_profile(muts, asm=None):
    cn = complementary_nucleotide
    mutations = defaultdict(float)
    N_skipped = 0

    raw_mutations = []
    for line in muts:
        if line.startswith("#"):
            continue
        if len(line) < 10:
            continue

        col_list = line.split()
        if len(col_list) < 4:
            continue

        # ID = col_list[2]

        # chromosome is expected to be one or two number or one letter
        chrom = col_list[0]  # VCF CHROM
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]
        # if len(chrom) == 2 and chrom[1] not in "0123456789":
        #     chrom = chrom[0]

        pos = int(col_list[1])  # VCF POS
        x = col_list[3]         # VCF REF
        y = col_list[4]         # VCF ALT

        # if multiple REF or ALT alleles are given, ignore mutation entry (could mean seq error, could mean deletion)
        if len(x) != 1:
            N_skipped += 1
            continue
        if len(y) != 1:
            N_skipped += 1
            continue

        raw_mutations.append((chrom, pos, x, y))
    # print("RAW", raw_mutations)

    MAX = 100000000
    if len(raw_mutations) > 0:
        if len(raw_mutations) > MAX:
            raw_mutations = raw_mutations[:MAX]

        contexts = get_context_batch(raw_mutations, asm)
        # print("CONTEXTS", contexts)

        if contexts is None:
            return None, None

        if len(contexts) == 0:
            return None, None

        for (chrom, pos, x, y) in raw_mutations:
            p5, p3 = contexts.get((chrom, pos), ("N", "N"))
            # print("RESULT: {} {}".format(p5, p3))

            if len(set([p5, x, y, p3]) - set(nucleotides)) > 0:
                # print(chrom, pos, p5, p3, x)
                # print("Skipping invalid nucleotides")
                N_skipped += 1
                continue

            if x in "CT":
                mutations[p5 + p3 + x + y] += 1.0
            else:
                # complementary mutation
                mutations[cn[p3] + cn[p5] + cn[x] + cn[y]] += 1.0

    N_loaded = int(sum(mutations.values()))
    processing_stats = {'loaded': N_loaded, 'skipped': N_skipped, 'format': 'VCF'}
    return mutations, processing_stats


def mp_ensembl_worker(raw_mutations_chunk):
    """ Helper function for ENSEMBL API call """

    cn = complementary_nucleotide
    verify_reference = {}

    contexts = {}  # result

    api_input = {'regions': []}
    for i, (assembly, chrom, pos, x, y) in enumerate(raw_mutations_chunk):
        # print("!!!", assembly, chrom, pos, x, y)
        api_input['regions'].append("{}:{}..{}".format(chrom, int(pos) - 1, int(pos) + 1))
        verify_reference[(chrom, int(pos))] = x
    api_input_json = json.dumps(api_input)

    asm = {
        38: "GRCh38",
        37: "GRCh37",
        36: "NCBI36",
        35: "NCBI35",
        34: "NCBI34"}

    server = "https://rest.ensembl.org"
    ext = "/sequence/region/human/?coord_system_version={}".format(asm.get(int(assembly), "GRCh38"))

    r = requests.post(server + ext, headers={
        "Content-Type": "application/json",
        "Accept": "application/json"}, data=api_input_json, timeout=2.000, verify=False)  # 2 second timeout
    if not r.ok:
        print("ENSEMBL Exception")
        return contexts
        # raise Exception("POST request does not work")
    decoded = r.json()

    if 'error' in decoded:
        return contexts

    for i, record in enumerate(decoded):
        if 'id' not in record:
            continue
        _, asm, chrom, begin, end, expand = record['id'].split(":")
        pos = int(begin) + 1
        if _ != "chromosome":
            continue

        seq = record['seq']
        nuc5 = seq[0]
        nuc = seq[1]
        nuc3 = seq[2]
        x = verify_reference.get((chrom, pos))
        if nuc != x:
            if cn[nuc] == x:
                nuc3 = cn[nuc5]
                nuc5 = cn[nuc3]
            else:
                # print("Could not match reference allele in ", raw_mutations_chunk[i], "in sequence", seq)
                nuc3 = 'N'
                nuc5 = 'N'
        # print("{} {}: {} [{} > ?] {}".format(chrom, pos, nuc5, x, nuc3))
        contexts[(chrom, pos)] = (nuc5, nuc3)
    return contexts


def get_context_ensembl(mutations, assembly):
    """
    Use Ensembl to get context of mutations
    ENSEMBL REST API http://rest.ensembl.org/documentation/info/sequence_region_post
    """
    raw_mutations = []
    for (chrom, pos, x, y) in mutations:
        raw_mutations.append((assembly, chrom, pos, x, y))

    L = len(raw_mutations)
    MAX_POST_SIZE = 50  # https://github.com/Ensembl/ensembl-rest/blob/78eebccb589798662c0ea758081ba2b2cb6acbca/ensembl_rest.conf.default     max_post_size = 50

    chunks = []
    for chunk in range(0, L, MAX_POST_SIZE):
        chunks.append(raw_mutations[chunk: chunk + MAX_POST_SIZE])

    try:
        p = multiprocessing.Pool(7)
        cc = p.map(mp_ensembl_worker, chunks, MAX_POST_SIZE)
    except:
        # raise
        cc = []
    finally:
        p.close()
        p.join()

    contexts = {}
    for c in cc:
        contexts.update(c)

    if len(contexts) != len(raw_mutations):
        # print("Not all mutations have context")
        pass
    return contexts


def read_cohort_mutations_from_tar(tar_fname, cohort):
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


def get_genome_sequence_twobit(twobit):
    cn = complementary_nucleotide
    for (chrom, pos, x, y) in mutations:
        start = int(pos) - 1  # zero-based numbering
        chrom = str(chrom)
        chromosome = chrom if chrom.startswith('chr') else 'chr' + chrom

        nuc5 = 'N'
        nuc3 = 'N'
        nuc = 'N'
        if chromosome in f:
            try:
                seq = f[chromosome][start - 1: start + 2]  # +/- 1 nucleotide
                nuc5, nuc, nuc3 = tuple(seq.upper())
            except:
                nuc = 'N'
                # print(chromosome, x, nuc5, nuc, nuc3)
            if nuc != 'N' and nuc != x:
                if cn[nuc] == x:
                    nuc3 = cn[nuc5]
                    nuc5 = cn[nuc3]
                else:
                    nuc3 = nuc5 = 'N'
        else:
            # print("NO CHROM", chromosome)
            pass
        contexts[(chrom, pos)] = (nuc5, nuc3)
    return contexts


def read_MAF_with_genomic_context_fname(fname, genome):
    with open(fname) as infile:
        read_MAF_with_genomic_context(infile, genome)


def read_MAF_with_genomic_context(infile, genome):
    twobit = tbr.TwoBitFile(genome)

    mutations = []
    reader = csv.reader(infile, delimiter='\t')
    # get names from column headers
    # MAF = namedtuple("MAF", map(str.lower, next(reader)))
    header = next(reader)
    # print(header)
    MAF = namedtuple("MAF", header)
    # print(MAF)
    for data in map(MAF._make, reader):
        if '_' in data.Protein_Change:
            continue
        if data.Protein_Change.endswith('fs'):
            continue
        if not data.Codon_Change.startswith("c."):
            continue
        c1, c2 = data.Codon_Change.split(")")[1].split(">")
        for offset in range(3):
            if c1[offset] != c2[offset]:
                break
        chrom = data.Chromosome
        start = int(data.Start_position)
        chromosome = chrom if chrom.startswith('chr') else 'chr' + chrom
        seq5_fwd = twobit[chromosome][start - offset - 2: start - offset + 3].upper()
        seq5_rev = "".join(
            [complementary_nucleotide[x] for x in reversed(
                twobit[chromosome][start - (2 - offset) - 2: start - (2 - offset) + 3].upper())])

        if seq5_fwd[1:4] == c1:
            seq5 = seq5_fwd
        elif seq5_rev[1:4] == c1:
            seq5 = seq5_rev
        else:
            print("Sequence missmatch of mutation with the genome, check if using the correct genome assembly:", data.Protein_Change)
            continue
        # data.Chromosome,
        # data.Start_position,
        # data.Strand,
        # data.Codon_Change,
        mutations.append((
            data.Hugo_Symbol,
            data.Protein_Change.split(".")[1],
            seq5
        ))

        # print(data.Codon_Change, data.Strand, offset, seq5)

        # print(data.cDNA_position, data.CDS_position, data.Protein_position)
        # print(data.Codon_Change)  #, data.CDS_position, data.Protein_position)
        # cDNA_position   Relative position of base pair in the cDNA sequence as a fraction. A "-" symbol is displayed as the numerator if the variant does not appear in cDNA
        # 54 - CDS_position   Relative position of base pair in coding sequence. A "-" symbol is displayed as the numerator if the variant does not appear in coding sequence
        # 55 - Protein_position
        # print(data.HGVSp)
        # Hugo_Symbol
        # Chromosome
        # Start_Position
        # End_Position
        # Strand
        # HGVSp
    return mutations


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


def fetch_genome(name):
    tbr.download.save_genome(name, destdir=None, mode='ftp')


def fetch_cohorts():
    urllib.request.urlretrieve("https://www.ncbi.nlm.nih.gov/mutagene/static/data/cohorts.tar.gz", filename='cohorts.tar.gz')

