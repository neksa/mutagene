# import csv
import numpy as np
from scipy import stats
from collections import defaultdict

# def write_motif_matches(outfile, motif_matches):
#     if len(motif_matches) > 0:
#         header = motif_matches[0].keys()
#         writer = csv.DictWriter(outfile, fieldnames=header, dialect='excel-tab')
#         writer.writeheader()
#         for data in motif_matches:
#             writer.writerow(data)


def write_multisample_decomposition(fname, samples_results, signature_ids, write_zeros=False):
    if isinstance(fname, (str, bytes)):
        o = open(fname, 'w')
    else:
        # fname is assumed to be a file descriptor
        o = fname

    fname.write("sample\tsignature\texposure\tmutations\n")

    for sample, results in samples_results.items():
        if type(results) is np.ndarray:
            h = results
        else:
            exposure_dict = {x['name']: x['score'] for x in results}
            exposure = [exposure_dict[name] for name in signature_ids]
            h = np.array(exposure)

            mutations_dict = {x['name']: x['mutations'] for x in results}
            mutations = [mutations_dict[name] for name in signature_ids]
            m = np.array(mutations, np.int)
        for i in range(h.shape[0]):
            if np.isclose(m[i], 0.0) and not write_zeros:
                continue
            fname.write("{}\t{}\t{:.4f}\t{}\n".format(sample, signature_ids[i], h[i], m[i]))

    # if filename then close descriptor
    if isinstance(fname, (str, bytes)):
        o.close()


def write_bootstrap_multisample_decomposition(fname, samples_results, signature_ids, write_zeros=False):
    if isinstance(fname, (str, bytes)):
        o = open(fname, 'w')
    else:
        # fname is assumed to be a file descriptor
        o = fname

    fname.write("sample\tsignature\texposure\tmutations\n")

    for sample, results in samples_results.items():
        if type(results) is np.ndarray:
            h = results
        else:
            exposure_dict = {x['name']: x['score'] for x in results}
            exposure = [exposure_dict[name] for name in signature_ids]
            h = np.array(exposure)

            mutations_dict = {x['name']: x['mutations'] for x in results}
            mutations = [mutations_dict[name] for name in signature_ids]
            m = np.array(mutations, np.int)
        for i in range(h.shape[0]):
            if np.isclose(m[i], 0.0) and not write_zeros:
                continue
            fname.write("{}\t{}\t{:.4f}\t{}\n".format(sample, signature_ids[i], h[i], m[i]))

    # if filename then close descriptor
    if isinstance(fname, (str, bytes)):
        o.close()


def write_decomposition(fname, results, signature_ids, sample_name="sample", write_zeros=False):
    if type(results) is np.ndarray:
        h = results
        m = np.round(h * 100)  # artificial number of mutations in case we only have a vector
    else:
        exposure_dict = {x['name']: x['score'] for x in results}
        exposure = [exposure_dict[name] for name in signature_ids]
        h = np.array(exposure)

        mutations_dict = {x['name']: x['mutations'] for x in results}
        mutations = [mutations_dict[name] for name in signature_ids]
        m = np.array(mutations, np.int)

    outf = open(fname, 'w') if isinstance(fname, (str, bytes)) else fname

    outf.write("sample\tsignature\texposure\tmutations\n")
    header = "{}\t{}\t{:.4f}\t{}\n"
    for i in range(h.shape[0]):
        if np.isclose(m[i], 0.0) and not write_zeros:
            continue
        outf.write(header.format(sample_name, signature_ids[i], h[i], m[i]))

    if isinstance(fname, (str, bytes)):
        outf.close()


def write_bootstrap_decomposition(fname, results, signature_ids, sample_name, write_zeros=False):
    exposures_lists = defaultdict(list)
    mutations_lists = defaultdict(list)
    for result in results:
        for x in result:
            exposures_lists[x['name']].append(x['score'])
            mutations_lists[x['name']].append(x['mutations'])

    # print(exposures_lists)
    # print(mutations_lists)

    h = np.array([np.nanmean(exposures_lists[name]) for name in signature_ids])
    m = np.array([np.nanmean(mutations_lists[name]) for name in signature_ids])  # , np.int)

    h_sem = np.array([stats.sem(exposures_lists[name]) for name in signature_ids])
    m_sem = np.array([stats.sem(mutations_lists[name]) for name in signature_ids])

    # FIXME: code duplication
    if isinstance(fname, (str, bytes)):
        o = open(fname, 'w')
    else:
        o = fname
    o.write("sample\tsignature\texposure\tmutations\texposure_CI_low\texposure_CI_high\tmutations_CI_low\tmutations_CI_high\n")
    for i in range(h.shape[0]):
        h_ci_low = h[i] - 1.96 * h_sem[i]
        h_ci_hi = h[i] + 1.96 * h_sem[i]
        m_ci_low = m[i] - 1.96 * m_sem[i]
        m_ci_hi = m[i] + 1.96 * m_sem[i]

        # if np.isclose(m[i], 0.0) and not write_zeros:
        #     continue
        if m_ci_low < 1.0:
            continue

        o.write("{}\t{}\t{:.4f}\t{:.0f}\t{:.4f}\t{:.4f}\t{:.0f}\t{:.0f}\n".format(sample_name, signature_ids[i], h[i], m[i], h_ci_low, h_ci_hi, m_ci_low, m_ci_hi))
    if isinstance(fname, (str, bytes)):
        o.close()


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
