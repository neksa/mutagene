import pathlib
import glob
import os
# import pprint
import numpy as np

from math import ceil
from os.path import isfile
from multiprocessing import Pool

from mutagene.io.profile import write_profile, read_profile_file
from mutagene.io.decomposition import write_decomposition, read_decomposition
from mutagene.signatures.identify import decompose_mutational_profile_counts
from mutagene.benchmark.deconstructsigs import deconstruct_sigs, deconstruct_sigs_custom


def gen_benchmark_2combinations(data_root, signature_names, W):
    name_to_idx = {}
    for i, s in enumerate(signature_names):
        name_to_idx[s] = i

    signature_names = list(name_to_idx.keys())
    # print(signature_names)
    # signature_ids = [str(name.split()[1]) for name in signature_names]

    # for r in [0.9, 0.8, 0.7, 0.6, 0.5]:
    for r in [0.7, 0.5]:
        for l in [0.1, 0.2]:
            # for n in [10, 50, 100, 500, 1000, 10000]:
            for n in [50, 500]:
                gen_sample_2combinations(data_root, signature_names, W, ratio=r, noise_level=l, n_mutations=n)


def gen_sample_2combinations(data_root, signature_ids, W, ratio=0.7, noise_level=0.05, n_mutations=10):
    assert 0.0 < ratio < 1.0
    assert 0.0 <= noise_level < ratio
    assert n_mutations > 0
    # np.set_printoptions(precision=4)
    N = W.shape[1]
    assert N >= 5

    n_replica = 10

    for i in range(N):
        for j in range(N):
            # if i == j:
            #     continue

            dirname = "{}/2comb/{}_{}_{}".format(data_root, N, signature_ids[i], signature_ids[j])
            pathlib.Path(dirname).mkdir(parents=True, exist_ok=True)
            print(dirname)

            for r in range(n_replica):
                fname = dirname + "/{:02d}_{:02d}_{:02d}_{:02d}".format(int(ratio * 10), int(noise_level * 100), int(n_mutations), r)
                # print(fname)

                profile_fname = fname + ".profile"
                info_fname = fname + ".info"

                h0 = np.zeros(N)
                h0[i] = ratio
                h0[j] = 1.0 - ratio
                h0 /= h0.sum()

                v0 = W.dot(h0)
                v0_counts = np.random.multinomial(ceil((1.0 - noise_level) * n_mutations), v0 / v0.sum())
                # print(v0_counts)
                v0_counts += np.random.multinomial(ceil(noise_level * n_mutations), [1.0 / v0.shape[0]] * v0.shape[0])  # uniform distribution
                # print(v0_counts)

                write_profile(profile_fname, v0_counts)
                # write_decomposition(info_fname, h0, signature_ids)
                write_decomposition(info_fname, h0, signature_ids, "synthetic", write_zeros=True)


def run_benchmark_2combinations(data_root, N, signature_ids, W, force=False):
    methods = ['MLE', 'MLEZ', 'AICc', 'BIC', 'AICcZ', 'BICZ']

    for fname in glob.glob("{}/2comb/{}_**/*.profile".format(data_root, N), recursive=True):
        print(fname)
        profile = read_profile_file(fname)

        for method in methods:
            info = "{}.{}.info".format(fname.split(".")[0], method)
            if isfile(info) and not force:
                continue

            _, _, results = decompose_mutational_profile_counts(
                profile,
                (W, signature_ids),
                method,
                debug=False,
                others_threshold=0.0)
            exposure_dict = {x['name']: x['score'] for x in results}
            exposure = [exposure_dict[name] for name in signature_ids]
            write_decomposition(info, np.array(exposure), signature_ids)


def run_benchmark_2combinations_deconstruct_sigs_helper(data):
    fname, ds_info, signature_ids, sigtype = data
    print(fname)
    results = deconstruct_sigs_custom(fname, signatures=int(sigtype))
    # print(results)
    exposure_dict = {x['name']: x['score'] for x in results}
    exposure = [exposure_dict[name] for name in signature_ids]
    write_decomposition(ds_info, np.array(exposure), signature_ids)


def run_benchmark_2combinations_deconstruct_sigs(data_root, N, signature_ids, W, force=False):
    def get_iterator():
        for fname in glob.glob("{}/2comb/{}_**/*.profile".format(data_root, N), recursive=True):
            sigtype = fname.split("/")[-2].split("_")[0]
            ds_info = fname.split(".")[0] + ".ds.info"
            if isfile(ds_info) and not force:
                continue
            yield (fname, ds_info, signature_ids, sigtype)

    with Pool(8) as p:
        p.map(run_benchmark_2combinations_deconstruct_sigs_helper, get_iterator(), 100)


def aggregate_benchmarks(data_root):
    methods = {
        "original": ".info",
        "mle": ".MLE.info",
        "mlez": ".MLEZ.info",
        "ds": ".ds.info",
        'aicc': '.AICc.info',
        'bic': '.BIC.info',
        'aiccz': '.AICcz.info',
        'bicz': '.BICz.info',
    }

    # generate full panel for all signatures in 30 x 30 signatures analysis
    with open("{}/2comb/res1-1.txt".format(data_root), 'w') as o:
        o.write("sigtype\tsig1\tsig2\tratio\tnoise\tnmut\treplica\tmethod\tsignature\tvalue\n")
        for fname in glob.glob("data/benchmark/2comb/**/*.profile", recursive=True):
            sigtype, sig1, sig2 = fname.split("/")[-2].split("_")
            ratio, noise, nmut, replica = fname.split("/")[-1].split(".")[0].split("_")

            for method in methods:
                info_fname = fname.split(".")[0] + methods[method]
                values, names = read_decomposition(info_fname)

                if values is None:
                    continue

                for value, signature in zip(values, names):
                    o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sigtype, sig1, sig2, ratio, noise, nmut, replica, method, signature, value))

    # only report the signature 2 value (as in DeconstructSigs benchmark)
    with open("{}/2comb/res2-2.txt".format(data_root), 'w') as o:
        o.write("sigtype\tsig1\tsig2\tratio\tnoise\tnmut\treplica\tmethod\tvalue\n")
        for fname in glob.glob("{}/2comb/**/*.profile".format(data_root), recursive=True):
            sigtype, sig1, sig2 = fname.split("/")[-2].split("_")
            ratio, noise, nmut, replica = fname.split("/")[-1].split(".")[0].split("_")

            for method in methods:
                # if method == 'ds' and sigtype != "30":
                #     continue

                info_fname = fname.split(".")[0] + methods[method]
                values, names = read_decomposition(info_fname)

                if values is None:
                    continue

                for value, signature in zip(values, names):
                    if sig2 == signature:
                        o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sigtype, sig1, sig2, ratio, noise, nmut, replica, method, value))


if __name__ == '__main__':
    pass
