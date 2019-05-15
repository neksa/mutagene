import multiprocessing
import json
import requests
from tqdm import tqdm

from mutagene.dna import complementary_nucleotide


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
        cc = p.map(mp_ensembl_worker, tqdm(chunks), MAX_POST_SIZE)
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
