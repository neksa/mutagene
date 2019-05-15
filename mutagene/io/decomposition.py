
import numpy as np


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
