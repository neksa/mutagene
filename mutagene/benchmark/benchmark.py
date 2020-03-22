#
import numpy as np

# from .io import read_profile
# from .io import format_profile
from mutagene.signatures.identify import decompose_mutational_profile_counts


def convert_to_list(name_to_idx, d):
    v = [0] * N
    for value in d:
        idx = name_to_idx.get(value['name'])
        if idx is None:
            # print(value['name'], 'not found in name_to_idx')
            continue
        v[idx] = value['score']
    return v


def get_scores(d):
    ll = 0
    frob = 0
    frob0 = 0
    js = 0
    kl = 0
    for value in d:
        if value['name'] == 'LogLik':
            ll = value['score']
        if value['name'] == 'Frobenius':
            frob = value['score']
        if value['name'] == 'FrobeniusZero':
            frob0 = value['score']
        if value['name'] == 'DivergenceJS':
            js = value['score']
        if value['name'] == 'DivergenceKL':
            kl = value['score']
    return ll, frob, frob0, js, kl


def benchmark_simulated(results_fname, signature_names, W):
    name_to_idx = {}
    for i, s in enumerate(signature_names):
        name_to_idx[s] = i

    W = np.array(W).T
    print(W.shape)
    N = W.shape[1]

    signature_names = list(name_to_idx.keys())
    signature_ids = [str(name.split()[1]) for name in signature_names]

    np.set_printoptions(precision=4)

    with open(results_fname, 'w') as report:
        report.write("TEST\tMUTATIONS\tMETH\tROUND\tSIGNATURE\tVALUE\tSE_E\tMSE_M\tMSE_E\tLL\tFROB\tFROB0\tJS\tKL\n")
        # for j in range(N):
        #     report.write("S{}\t".format(j + 1))
        # report.write("MSE\n")

        for i in range(N):

            for N_mutations in [10, 50, 100, 500, 1000, 10000]:

                for iteration in range(10):
                    h0 = np.zeros(N)
                    h0[i] = int(0.8 * N_mutations)
                    for j in range(N):
                        if i == j:
                            continue

                    # 20% uniform noise
                    # for k in range(int(0.2 * N_mutations)):
                    #     h0[random.randrange(N)] += 1
                    h0 += np.random.multinomial(int(0.2 * N_mutations), [1.0 / N] * N)  # uniform distribution
                    h0_counts = h0.copy()

                    h0 /= h0.sum()
                    v0 = W.dot(h0)
                    v0_counts = np.random.multinomial(N_mutations, v0 / v0.sum())  # np.ceil(v0 * N_mutations)

                    # oname = prefix + "profile"
                    # # print("ONAME:", oname)
                    # with open(oname, 'w') as o:
                    #     query_formatted = format_profile(h0.tolist())
                    #     o.write(query_formatted)

                    # oname = prefix + "counts"
                    # with open(oname, 'w') as o:
                    #     query_formatted = format_profile(h0.tolist(), counts=True)
                    #     o.write(query_formatted)

                    THRESHOLD = 0.0
                    # print(h0_counts)
                    # print(v0_counts)
                    _, _, exposure = decompose_mutational_profile_counts(v0_counts, W, 'MLE', others_threshold=THRESHOLD)

                    h1 = np.array(convert_to_list(name_to_idx, exposure))
                    # h1 += min(1.0 - h1.sum(), 0)
                    h1_counts = np.random.multinomial(N_mutations, h1 / h1.sum())
                    ll, frob, frob0, js, kl = get_scores(exposure)

                    # print(np.round(h0, 3))
                    # print(np.round(h1, 3))
                    v1 = W.dot(h1)

                    # convert exposure to mutational burden
                    # v1_counts = np.ceil(v1 * N_mutations)

                    # convert exposure to mutational burden
                    v1_counts = np.random.multinomial(N_mutations, v1 / v1.sum())

                    MSE_M = np.mean((v0_counts - v1_counts)**2)
                    MSE_E = np.mean((h0 - h1)**2)

                    for j in range(N):
                        report.write("{}_{}_{}\t{}\t{}\t{}\t".format(signature_ids[i], 80, 20, N_mutations, "MLE", iteration + 1))
                        report.write("{}\t{}\t".format(signature_names[j], int(h1_counts[j])))
                        SE_E = (h0[j] - h1[j])**2
                        report.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(SE_E, MSE_M, MSE_E, ll, frob, frob0, js, kl))

                    print(signature_ids[i], N_mutations, iteration, round(MSE_M, 4), round(MSE_E, 4), round(ll, 4), round(frob, 4), round(frob0, 4), round(js, 4), round(kl, 4))


def benchmark_2combinations(results_fname, signature_names, W):
    name_to_idx = {}
    for i, s in enumerate(signature_names):
        name_to_idx[s] = i

    W = np.array(W).T
    N = W.shape[1]

    signature_names = list(name_to_idx.keys())
    signature_ids = [str(name.split()[1]) for name in signature_names]

    np.set_printoptions(precision=4)

    noise_level = 0.05
    ratio = 0.7
    n_sample = 20

    with open(results_fname, 'w') as report:
        report.write("ratio\tn_sample\tnoise_level\tsignatures\tSIGNATURE\tVALUE\tSE_E\tMSE_M\tMSE_E\tLL\tFROB\tFROB0\tJS\tKL\n")
        # for j in range(N):
        #     report.write("S{}\t".format(j + 1))
        # report.write("MSE\n")

        for i in range(N):

            # for N_mutations in [10, 50, 100, 500, 1000, 10000]:
            for N_mutations in [50, 500]:

                for iteration in range(10):
                    h0 = np.zeros(N)
                    h0[i] = int(0.8 * N_mutations)
                    for j in range(N):
                        if i == j:
                            continue

                    # 20% uniform noise
                    # for k in range(int(0.2 * N_mutations)):
                    #     h0[random.randrange(N)] += 1
                    h0 += np.random.multinomial(int(0.2 * N_mutations), [1.0 / N] * N)  # uniform distribution
                    h0_counts = h0.copy()

                    h0 /= h0.sum()
                    v0 = W.dot(h0)
                    v0_counts = np.random.multinomial(N_mutations, v0 / v0.sum())  # np.ceil(v0 * N_mutations)

                    THRESHOLD = 0.0
                    # print(h0_counts)
                    # print(v0_counts)
                    _, _, exposure = decompose_mutational_profile_counts(v0_counts, W, 'MLE', others_threshold=THRESHOLD)

                    h1 = np.array(convert_to_list(name_to_idx, exposure))
                    # h1 += min(1.0 - h1.sum(), 0)
                    h1_counts = np.random.multinomial(N_mutations, h1 / h1.sum())
                    ll, frob, frob0, js, kl = get_scores(exposure)

                    # print(np.round(h0, 3))
                    # print(np.round(h1, 3))
                    v1 = W.dot(h1)

                    # convert exposure to mutational burden
                    # v1_counts = np.ceil(v1 * N_mutations)

                    # convert exposure to mutational burden
                    v1_counts = np.random.multinomial(N_mutations, v1 / v1.sum())

                    MSE_M = np.mean((v0_counts - v1_counts)**2)
                    MSE_E = np.mean((h0 - h1)**2)

                    for j in range(N):
                        report.write("{}_{}_{}\t{}\t{}\t{}\t".format(signature_ids[i], 80, 20, N_mutations, "MLE", iteration + 1))
                        report.write("{}\t{}\t".format(signature_names[j], int(h1_counts[j])))
                        SE_E = (h0[j] - h1[j])**2
                        report.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(SE_E, MSE_M, MSE_E, ll, frob, frob0, js, kl))

                    print(signature_ids[i], N_mutations, iteration, round(MSE_M, 4), round(MSE_E, 4), round(ll, 4), round(frob, 4), round(frob0, 4), round(js, 4), round(kl, 4))
