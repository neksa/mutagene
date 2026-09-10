import glob
import logging
import os
import pathlib
import random
import uuid
from multiprocessing import Pool

import numpy as np
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    mean_squared_error,
    precision_score,
    recall_score,
)

from mutagene.benchmark.deconstructsigs import deconstruct_sigs_custom
from mutagene.benchmark.generate_benchmark import *
from mutagene.io.profile import read_profile_file, read_signatures, write_profile
from mutagene.signatures.identify import NegLogLik

logger = logging.getLogger(__name__)

# from mutagene.identify import decompose_mutational_profile_counts


def multiple_benchmark_helper(task):
    """Generate one synthetic sample and decompose it.

    Takes its directory and signature set as arguments. They used to be
    hardcoded here, so --root and --signatures were accepted by the command
    line and then ignored.
    """
    j, dirname, signature_sets, run_deconstruct_sigs = task

    for name in signature_sets:
        W, signature_names = read_signatures(name)
        N = W.shape[1]

        # How many signatures to mix. Bounded by how many the set actually has:
        # this used the set's name, which is a string like "COSMICv3".
        r = random.randrange(2, min(N + 1, 15))

        # print(np.random.choice(N, r), .05 + np.random.dirichlet(np.ones(r), 1))
        while True:
            h0 = np.zeros(N)
            h0[np.random.choice(N, r)] = 0.05 + np.random.dirichlet(np.ones(r), 1)
            if np.greater(h0, 0.05).sum() == r:
                break
        h0 /= h0.sum()
        v0 = W.dot(h0)
        # print(h0)
        n_mutations = random.randrange(10, 50)
        v0_counts = np.random.multinomial(n_mutations, v0 / v0.sum())
        # print(v0_counts)

        random_name = str(uuid.uuid4())[:4]
        fname = os.path.join(dirname, f"{name}_{r}_{n_mutations}_{random_name}")
        print(fname)
        profile_fname = fname + ".profile"
        info_fname = fname + ".info"
        mle_info = fname + ".MLE.info"
        mlez_info = fname + ".MLEZ.info"
        ds_info = fname + ".ds.info"

        write_profile(profile_fname, v0_counts)
        write_decomposition(info_fname, {"synthetic": h0}, signature_names)

        # deconstructSigs is an R package, so this is only attempted when the
        # caller asked for it and Rscript is actually installed.
        if run_deconstruct_sigs:
            results = deconstruct_sigs_custom(profile_fname, signatures=name)
            write_decomposition(ds_info, {"synthetic": results}, signature_names)
        profile = read_profile_file(profile_fname)
        for method, method_fname in [("MLE", mle_info), ("MLEZ", mlez_info)]:
            _, _, results = decompose_mutational_profile_counts(
                profile, (W, signature_names), method, others_threshold=0.0
            )
            write_decomposition(method_fname, {"synthetic": results}, signature_names)


def multiple_benchmark(
    data_root="data/benchmark",
    signature_sets=(30,),
    replicates=100,
    processes=10,
    run_deconstruct_sigs=False,
):
    """Generate synthetic multi-signature samples under data_root/multiple."""
    dirname = os.path.join(data_root, "multiple")
    pathlib.Path(dirname).mkdir(parents=True, exist_ok=True)
    random.seed(13425)

    tasks = [(j, dirname, list(signature_sets), run_deconstruct_sigs) for j in range(replicates)]
    if processes == 1:
        # A pool of one is all overhead, and it makes failures much harder to
        # read because the traceback comes back through the parent.
        for task in tasks:
            multiple_benchmark_helper(task)
    else:
        with Pool(processes) as p:
            p.map(multiple_benchmark_helper, tasks)


def multiple_benchmark_run_helper(data):
    fname, signature_ids, W, force = data
    # methods = ['MLE', 'MLEZ', 'AICc', 'BIC', 'AICcZ', 'BICZ']
    methods = ["AICc", "AICcZ"]

    # print(fname)
    profile = read_profile_file(fname)

    for method in methods:
        info = "{}.{}.info".format(fname.split(".")[0], method)
        if isfile(info) and not force:
            continue

        print(info)

        _, _, results = decompose_mutational_profile_counts(
            profile, (W, signature_ids), method, others_threshold=0.0
        )
        exposure_dict = {x["name"]: x["score"] for x in results}
        exposure = [exposure_dict[name] for name in signature_ids]
        write_decomposition(info, {"synthetic": np.array(exposure)}, signature_ids)


def multiple_benchmark_run(
    N, signature_ids, W, force=False, data_root="data/benchmark", processes=10
):
    """Decompose every generated profile for signature set N."""
    pattern = os.path.join(data_root, "multiple", f"{N}_*.profile")
    tasks = [(fname, signature_ids, W, force) for fname in sorted(glob.glob(pattern))]
    if not tasks:
        logger.warning(
            f"No generated profiles found in {pattern}. Run the matching *_gen mode first"
        )
        return

    random.seed(13425)
    if processes == 1:
        for task in tasks:
            multiple_benchmark_run_helper(task)
    else:
        with Pool(processes) as p:
            p.map(multiple_benchmark_run_helper, tasks)


def aggregate_multiple_benchmarks(data_root="data/benchmark"):
    methods = {
        "mle": ".MLE.info",
        "mlez": ".MLEZ.info",
        "ds": ".ds.info",
        "aicc": ".AICc.info",
        "bic": ".BIC.info",
        "aiccz": ".AICcz.info",
        "bicz": ".BICz.info",
    }

    # signatures_thresholds = {
    #     5: 0.06,
    #     10: 0.03,
    #     30: 0.01,
    # }

    # signatures_thresholds = {
    #     5: 0.0001,
    #     10: 0.0001,
    #     30: 0.0001,
    # }

    # only report the signature 2 value (as in DeconstructSigs benchmark)
    multiple_dir = os.path.join(data_root, "multiple")
    with open(os.path.join(multiple_dir, "res1.txt"), "w") as o:
        o.write(
            "file_id\tsigtype\tnsig\tnmut\tmethod\tSRMSE\tPRMSE\tSTRMSE\tLLIK\tLLIK0\tTLLIK\tTLLIK0\tprecision\trecall\taccuracy\tf1\n"
        )
        for fname in sorted(glob.glob(os.path.join(multiple_dir, "*.profile"))):
            file_id = fname.split("/")[-1].split(".")[0]
            sigtype, r, nmut, replica = fname.split("/")[-1].split(".")[0].split("_")
            sigtype = int(sigtype)

            if sigtype != 30:
                continue

            W, signature_names = read_signatures(sigtype)

            info_fname = fname.split(".")[0] + ".info"
            orig_profile = read_profile_file(fname)
            h0, names = read_decomposition(info_fname)

            # threshold = 0.06
            threshold = 0.06

            # threshold = 1.0 / np.sqrt(int(nmut)) if method != "ds" else 0.06
            h0_threshold = np.where(h0 > threshold, h0, 0.0)  # zero below threshold
            h0_binary = np.array(h0_threshold) > 0.0  # true / false for threshold
            nsig = np.count_nonzero(h0_binary)

            if nsig < int(r):
                print("LESS", sigtype, nsig, r)

            if nsig > int(r):
                print("MORE", sigtype, nsig, r)

            if nsig <= 1:
                continue
            if nsig > 10:
                continue

            for method in methods:
                method_fname = fname.split(".")[0] + methods[method]
                values, names = read_decomposition(method_fname)

                # print(method_fname)

                if values is None:
                    continue

                h = np.array(values)
                if h.sum() == 0:
                    continue

                h_threshold = np.where(h > threshold, h, 0.0)  # zero below threshold

                reconstructed_profile = W.dot(h)
                # print(h)
                # print(reconstructed_profile)

                PRMSE = np.sqrt(
                    mean_squared_error(
                        np.array(orig_profile) / np.array(orig_profile).sum(),
                        np.array(reconstructed_profile) / np.array(reconstructed_profile).sum(),
                    )
                )
                SRMSE = np.sqrt(mean_squared_error(h0, h))
                STRMSE = np.sqrt(mean_squared_error(h0_threshold, h_threshold))
                LLIK0 = -NegLogLik(h0, W, orig_profile)
                TLLIK0 = -NegLogLik(h0_threshold, W, orig_profile)
                LLIK = -NegLogLik(h, W, orig_profile)
                TLLIK = -NegLogLik(h_threshold, W, orig_profile)

                # print(h0.sum())
                # print(h.sum())

                h_binary = np.array(h_threshold) > 0.0  # true / false for threshold
                precision = precision_score(h0_binary, h_binary)
                recall = recall_score(h0_binary, h_binary)
                accuracy = accuracy_score(h0_binary, h_binary)
                f1 = f1_score(h0_binary, h_binary)

                o.write(
                    f"{file_id}\t{sigtype}\t{nsig}\t{nmut}\t{method}\t{SRMSE}\t{PRMSE}\t{STRMSE}\t{LLIK}\t{LLIK0}\t{TLLIK}\t{TLLIK0}\t{precision}\t{recall}\t{accuracy}\t{f1}\n"
                )
