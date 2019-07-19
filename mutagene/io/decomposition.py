# import csv
import numpy as np

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


def write_decomposition(fname, results, signature_ids, sample_name, write_zeros=False):
    if type(results) is np.ndarray:
        h = results
    else:
        exposure_dict = {x['name']: x['score'] for x in results}
        exposure = [exposure_dict[name] for name in signature_ids]
        h = np.array(exposure)

        mutations_dict = {x['name']: x['mutations'] for x in results}
        mutations = [mutations_dict[name] for name in signature_ids]
        m = np.array(mutations, np.int)

    # FIXME: code duplication
    if isinstance(fname, (str, bytes)):
        with open(fname, 'w') as o:
            fname.write("sample\tsignature\texposure\tmutations\n")
            for i in range(h.shape[0]):
                if np.isclose(m[i], 0.0) and not write_zeros:
                    continue
                o.write("{}\t{}\t{:.4f}\t{}\n".format(sample_name, signature_ids[i], h[i], m[i]))
    else:
        fname.write("sample\tsignature\texposure\tmutations\n")
        for i in range(h.shape[0]):
            if np.isclose(m[i], 0.0) and not write_zeros:
                continue
            fname.write("{}\t{}\t{:.4f}\t{}\n".format(sample_name, signature_ids[i], h[i], m[i]))


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
