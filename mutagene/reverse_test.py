import pathlib
import glob
import pprint
import random
import numpy as np
import uuid

# from math import ceil
# from os.path import isfile
# from multiprocessing import Pool

from .io import read_profile, format_profile, read_signatures
from .io import write_profile
# from .identify import decompose_mutational_profile_counts
from .deconstructsigs import deconstruct_sigs_custom
from .generate_benchmark import *


def reverse_test():
    dirname = "data/benchmark/reverse/"
    pathlib.Path(dirname).mkdir(parents=True, exist_ok=True)

    random.seed(13425)

    j = 0
    while j < 100:
        for i in [5, 10, 30]:
            W, signature_names = read_signatures(i)
            N = W.shape[1]

            random_name = str(uuid.uuid4())[:8]
            fname = dirname + "/{}".format(random_name)
            profile_fname = fname + ".profile"
            info_fname = fname + ".info"
            # mlez_info_f = fname + ".MLEZ-F.info"
            # mlez_info_r = fname + ".MLEZ-R.info"
            ds_info_f = fname + ".ds-F.info"
            ds_info_r = fname + ".ds-R.info"

            r = random.randrange(2, i // 2 + 1)
            h0 = np.zeros(N)
            h0[np.random.choice(N, r)] = 0.05 + np.round(np.random.dirichlet(np.ones(r), 1))
            h0 /= h0.sum()
            v0 = W.dot(h0)
            # print(h0)
            n_mutations = random.randrange(50, 1000)
            v0_counts = np.random.multinomial(n_mutations, v0 / v0.sum())
            # print(v0_counts)

            write_profile(profile_fname, v0_counts)
            write_decomposition(info_fname, h0, signature_names)

            ##################################################
            for direction, ds_fname in [(False, ds_info_f), (True, ds_info_r)]:
                results = deconstruct_sigs_custom(profile_fname, signatures=i, reverse=direction)
                exposure_dict = {x['name']: x['score'] for x in results}
                exposure = [exposure_dict[name] for name in signature_names]
                write_decomposition(ds_fname, np.array(exposure), signature_names)
        j += 1
