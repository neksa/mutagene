import json
from subprocess import Popen, PIPE, TimeoutExpired


def deconstruct_sigs(profile_fname, sample='sample'):
    script1 = """
library(jsonlite)
library(deconstructSigs)
"""

    # # if deconstruct_sigs.proc is None:
    # if 'proc' not in deconstruct_sigs.__dict__ or deconstruct_sigs.proc is None:
    #     proc = deconstruct_sigs.proc = Popen(["Rscript", "-"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    #     proc.stdin.write(script1.encode("utf-8"))
    #     out = proc.stdout.read()
    # else:
    #     proc = deconstruct_sigs.proc


    script2 = """s <- t(read.table('{}', sep="\t", header=FALSE, row.names=1))
row.names(s) <- '{}'
s <- as.data.frame(s)
w <- whichSignatures(s / sum(s), signatures.ref=signatures.cosmic, signature.cutoff=0.00)
toJSON(w)
"""
    script2 = (script1 + script2).format(profile_fname, sample).encode("utf-8")

    proc = Popen(["Rscript", "-"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    out, err = proc.communicate(script2, timeout=10)
    # exitcode = proc.returncode
    # print(exitcode, out, err)
    json_string = out.decode("utf-8")
    w = json.loads(json_string)
    # import pprint
    # pprint.pprint(w)
    result = []
    for k, v in w['weights'][0].items():
        if k.startswith('_row'):
            continue
        # if float(v) == 0.0:
        #     continue
        name = k.replace('Signature.', '')
        result.append({
            'name': name,
            'score': v})
    return result


def deconstruct_sigs_custom(profile_fname, sample='sample', signatures=30, cutoff=0.00):
    script1 = """
library(jsonlite)
library(deconstructSigs)
"""

    script1 += "W <- data.frame()"

    sig_map = {
        5: "A",
        10: "B",
        30: "C"
    }

    for i in range(signatures):
        script1 += """
W <- rbind(W, t(read.table('/Users/agoncear/projects/pymutagene/data/signatures/{}_{}.profile', sep="\t", header=FALSE, row.names=1)))
""".format(sig_map[signatures], i + 1)

    script1 += "row.names(W) <- c(" + ",".join(["'Signature." + str(i + 1) + "'" for i in range(signatures)]) + ")\n"
    script1 += "W <- as.data.frame(W)\n"

    script2 = """s <- t(read.table('{}', sep="\t", header=FALSE, row.names=1))
row.names(s) <- '{}'
s <- as.data.frame(s)
w <- whichSignatures(s / sum(s), signatures.ref=W, signature.cutoff={})
toJSON(w)
"""
    script2 = (script1 + script2).format(profile_fname, sample, cutoff)
    # print(script2)
    script2 = script2.encode("utf-8")

    proc = Popen(["Rscript", "-"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    out, err = proc.communicate(script2, timeout=10)
    # exitcode = proc.returncode
    # print(exitcode, out, err)
    json_string = out.decode("utf-8")
    w = json.loads(json_string)
    # import pprint
    # pprint.pprint(w)
    result = []
    for k, v in w['weights'][0].items():
        if k.startswith('_row'):
            continue
        # if float(v) == 0.0:
        #     continue
        name = k.replace('Signature.', '')
        result.append({
            'name': name,
            'score': v})
    return result
