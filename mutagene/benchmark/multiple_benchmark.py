import glob
import random
import uuid
import numpy as np

from multiprocessing import Pool
from sklearn.metrics import (
    recall_score, precision_score, accuracy_score, f1_score, mean_squared_error)

from mutagene.io.profile import read_profile_file, write_profile, read_signatures
from mutagene.signatures.identify import NegLogLik
from mutagene.benchmark.deconstructsigs import deconstruct_sigs_custom
from mutagene.benchmark.generate_benchmark import *
# from mutagene.identify import decompose_mutational_profile_counts


def multiple_benchmark_helper(j):
    dirname = "data/benchmark/multiple"

    # for i in [5, 10, 30]:
    for i in [30, ]:
        W, signature_names = read_signatures(i)
        N = W.shape[1]

        # r = random.randrange(2, i // 3 + 2)
        r = random.randrange(2, min(i + 1, 15))

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
        fname = dirname + "/{:02d}_{}_{}_{}".format(i, r, n_mutations, random_name)
        print(fname)
        profile_fname = fname + ".profile"
        info_fname = fname + ".info"
        mle_info = fname + ".MLE.info"
        mlez_info = fname + ".MLEZ.info"
        ds_info = fname + ".ds.info"

        write_profile(profile_fname, v0_counts)
        write_decomposition(info_fname, h0, signature_names)

        ##################################################
        results = deconstruct_sigs_custom(profile_fname, signatures=i)
        write_decomposition(ds_info, results, signature_names)
        ##################################################
        profile = read_profile_file(profile_fname)
        for method, method_fname in [("MLE", mle_info), ("MLEZ", mlez_info)]:
            _, _, results = decompose_mutational_profile_counts(
                profile,
                (W, signature_names),
                method,
                debug=False,
                others_threshold=0.0)
            write_decomposition(method_fname, results, signature_names)


def multiple_benchmark():
    # pathlib.Path(dirname).mkdir(parents=True, exist_ok=True)
    random.seed(13425)

    with Pool(10) as p:
        p.map(multiple_benchmark_helper, range(100))


def multiple_benchmark_run_helper(data):
    fname, signature_ids, W, force = data
    # methods = ['MLE', 'MLEZ', 'AICc', 'BIC', 'AICcZ', 'BICZ']
    methods = ['AICc', 'AICcZ']

    # print(fname)
    profile = read_profile_file(fname)

    for method in methods:
        info = "{}.{}.info".format(fname.split(".")[0], method)
        if isfile(info) and not force:
            continue

        print(info)

        _, _, results = decompose_mutational_profile_counts(
            profile,
            (W, signature_ids),
            method,
            debug=False,
            others_threshold=0.0)
        exposure_dict = {x['name']: x['score'] for x in results}
        exposure = [exposure_dict[name] for name in signature_ids]
        write_decomposition(info, np.array(exposure), signature_ids)


def multiple_benchmark_run(N, signature_ids, W, force=False):
    def get_iterator():
        for fname in glob.glob("data/benchmark/multiple/{:02d}_*.profile".format(N), recursive=True):
                yield (fname, signature_ids, W, force)

    random.seed(13425)
    with Pool(10) as p:
        p.map(multiple_benchmark_run_helper, get_iterator(), 100)


def aggregate_multiple_benchmarks():
    methods = {
        "mle": ".MLE.info",
        "mlez": ".MLEZ.info",
        "ds": ".ds.info",
        'aicc': '.AICc.info',
        'bic': '.BIC.info',
        'aiccz': '.AICcz.info',
        'bicz': '.BICz.info',
    }

    # signatures_thresholds = {
    #     5: 0.06,
    #     10: 0.03,
    #     30: 0.01,
    # }

    signatures_thresholds = {
        5: 0.06,
        10: 0.06,
        30: 0.06,
    }

    # signatures_thresholds = {
    #     5: 0.0001,
    #     10: 0.0001,
    #     30: 0.0001,
    # }

    # only report the signature 2 value (as in DeconstructSigs benchmark)
    with open("data/benchmark/multiple/res1.txt", 'w') as o:
        o.write("file_id\tsigtype\tnsig\tnmut\tmethod\tSRMSE\tPRMSE\tSTRMSE\tLLIK\tLLIK0\tTLLIK\tTLLIK0\tprecision\trecall\taccuracy\tf1\n")
        for fname in glob.glob("data/benchmark/multiple/*.profile", recursive=True):
            file_id = fname.split("/")[-1].split(".")[0]
            sigtype, r, nmut, replica = fname.split("/")[-1].split(".")[0].split("_")
            sigtype = int(sigtype)

            if sigtype != 30:
                continue

            W, signature_names = read_signatures(sigtype)

            info_fname = fname.split(".")[0] + '.info'
            orig_profile = read_profile_file(fname)
            h0, names = read_decomposition(info_fname)

            # threshold = 0.06
            threshold = 0.06

            # threshold = 1.0 / np.sqrt(int(nmut)) if method != "ds" else 0.06
            h0_threshold = np.where(h0 > threshold, h0, 0.0)  # zero below threshold
            h0_binary = np.array(h0_threshold) > 0.0   # true / false for threshold
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

                PRMSE = np.sqrt(mean_squared_error(
                    np.array(orig_profile) / np.array(orig_profile).sum(),
                    np.array(reconstructed_profile) / np.array(reconstructed_profile).sum()))
                SRMSE = np.sqrt(mean_squared_error(h0, h))
                STRMSE = np.sqrt(mean_squared_error(h0_threshold, h_threshold))
                LLIK0 = - NegLogLik(h0, W, orig_profile)
                TLLIK0 = - NegLogLik(h0_threshold, W, orig_profile)
                LLIK = - NegLogLik(h, W, orig_profile)
                TLLIK = - NegLogLik(h_threshold, W, orig_profile)

                # print(h0.sum())
                # print(h.sum())

                h_binary = np.array(h_threshold) > 0.0  # true / false for threshold
                precision = precision_score(h0_binary, h_binary)
                recall = recall_score(h0_binary, h_binary)
                accuracy = accuracy_score(h0_binary, h_binary)
                f1 = f1_score(h0_binary, h_binary)

                o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    file_id, sigtype, nsig, nmut, method, SRMSE, PRMSE, STRMSE, LLIK, LLIK0, TLLIK, TLLIK0, precision, recall, accuracy, f1))
