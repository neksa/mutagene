# import csv
import numpy as np
from scipy import stats
from collections import defaultdict
import pandas as pd


def _get_stats(results):
    """
        Convert dictionary decomposition results into a pandas dataframe compatible with bootstrap stats results
    """
    if type(results) is np.ndarray:
        h = results
    else:
        exposure_dict = {x['name']: x['score'] for x in results if x['mutations'] != ''}
        mutations_dict = {x['name']: x['mutations'] for x in results if x['mutations'] != ''}
        signatures = list(exposure_dict.keys())

        exposure = [exposure_dict[name] for name in signatures]
        mutations = [mutations_dict[name] for name in signatures]

        h = np.array(exposure)
        m = np.array(mutations, np.int)

    df = pd.DataFrame(
        dict(
            zip(
                "signature\texposure\tmutations".split(),
                [signatures, h, m])))
    return df.sort_values(by=["mutations", "exposure"], ascending=False)


def _get_bootstrap_stats_percentile(bootstrap_results, level):
    """
        level is converted to percentile: level 90 means an interval 9% - 95%

        Quantile-based bootstrap intervals
        Returns pandas data frame - each row is a signature, zeros not removed
    """
    ci_high = (100.0 + level) / 2.0
    ci_low = (100.0 - level) / 2.0

    exposures_lists = defaultdict(list)
    mutations_lists = defaultdict(list)
    signatures = set()

    for results in bootstrap_results:
        for x in results:
            if x['mutations'] != '' and x['score'] != '':
                exposures_lists[x['name']].append(x['score'])
                mutations_lists[x['name']].append(x['mutations'])
                signatures |= set([x['name'], ])
    signatures = list(signatures)

    # h: exposures, float [0..1]
    exposures = np.array([exposures_lists[name] for name in signatures]).T

    h = np.percentile(exposures, 50.0, axis=0)  # median
    h_ci_low = np.percentile(exposures, ci_low, axis=0)
    h_ci_high = np.percentile(exposures, ci_high, axis=0)

    # m: mutations, integer [0, 1, 2...]
    mutations = np.array([mutations_lists[name] for name in signatures]).T
    m = np.round(np.percentile(mutations, 50.0, axis=0))
    m_ci_low = np.round(np.percentile(mutations, ci_low, axis=0))
    m_ci_high = np.round(np.percentile(mutations, ci_high, axis=0))

    df = pd.DataFrame(
        dict(
            zip(
                "signature\texposure\tmutations\texposure_low\texposure_high\tmutations_low\tmutations_high".split(),
                [signatures, h, m, h_ci_low, h_ci_high, m_ci_low, m_ci_high])))
    return df.sort_values(by=["mutations", "exposure"], ascending=False)


def _get_bootstrap_stats_t(sample_results, bootstrap_results, n, level):
    """
        Calculate approximate bias-corrected t-based confidence intervals for a list of samples using Manly (2007) method.

        n is the number of mutations

        t-based confidence intervals for given confidence level will be based on t distribution with n-1 degrees of freedom

        bias-corrected means are calculated as 2 * theta - mean(bootstrap_theta), where:
            - theta is exposure calculated without the bootstrap
            - mean(bootstrap_theta) is expected value of exposure in bootstrap samples

        Returns pandas data frame - each row is a signature, zeros not removed
    """
    ci_level = (1.0 + level / 100.0) / 2.0

    # ppf: Percent point function (inverse of cdf — percentiles) of t-distribution with (n-1) d.f.
    t = stats.t.ppf(ci_level, n - 1)

    exposure = defaultdict(float)
    mutations = defaultdict(float)
    exposures_lists = defaultdict(list)
    mutations_lists = defaultdict(list)
    signatures = set()

    for x in sample_results:
        if x['mutations'] != '' and x['score'] != '':
            exposure[x['name']] = x['score']
            mutations[x['name']] = x['mutations']

    for results in bootstrap_results:
        for x in results:
            if x['mutations'] != '' and x['score'] != '':
                exposures_lists[x['name']].append(x['score'])
                mutations_lists[x['name']].append(x['mutations'])
                signatures |= set([x['name'], ])
    signatures = list(signatures)

    # h: exposures, float [0..1]
    h = np.array([2 * exposure[name] - np.nanmean(exposures_lists[name]) for name in signatures])
    h_sem = np.array([stats.sem(exposures_lists[name]) for name in signatures])
    h_ci_low = np.clip(h - t * h_sem, a_min=0.0, a_max=None)
    h_ci_high = np.clip(h + t * h_sem, a_min=None, a_max=1.0)
    h = np.clip(h, a_min=0.0, a_max=None)

    # m: mutations, integer [0, 1, 2...]
    m = np.array([2 * mutations[name] - np.nanmean(mutations_lists[name]) for name in signatures])
    m_sem = np.array([stats.sem(mutations_lists[name]) for name in signatures])
    m_ci_low = np.clip(m - t * m_sem, a_min=0, a_max=None)
    m_ci_high = m + t * m_sem
    m = np.clip(m, a_min=0, a_max=None)

    df = pd.DataFrame(
        dict(
            zip(
                "signature\texposure\tmutations\texposure_low\texposure_high\tmutations_low\tmutations_high".split(),
                [signatures, h, m, h_ci_low, h_ci_high, m_ci_low, m_ci_high]
            )))
    return df.sort_values(by="exposure", ascending=False)


def write_decomposition(
    fname, samples_results, signature_ids, mutations_threshold=0,
    bootstrap_method=None, profile=None,
    bootstrap_results=None, bootstrap_level=None
):
    """
        Process sample results for the input data and (optionally) bootstrapped data
        and save a resulting table with or without the optional bootstrap columns
    """
    bootstrap = bootstrap_results is not None and bootstrap_level is not None and profile is not None

    # print("Bootstrap", bootstrap)
    dfs = []
    for sample in samples_results.keys():
        if bootstrap:
            if bootstrap_method == 't':
                # t-distribution-based confidence intervals
                p = profile[sample]
                n = int(np.sum(p))
                assert n > 1
                df_sample = _get_bootstrap_stats_t(samples_results[sample], bootstrap_results[sample], n, bootstrap_level)
            elif bootstrap_method == 'p':
                # percentile-based confidence intervals
                df_sample = _get_bootstrap_stats_percentile(bootstrap_results[sample], bootstrap_level)
            else:
                raise ValueError("Incorrect bootstrap_method value, only 't' or 'p' are recognized")
        else:
            df_sample = _get_stats(samples_results[sample])
        df_sample.insert(0, 'sample', sample)
        dfs.append(df_sample)

    # join results for multiple samples
    if len(dfs) == 0:
        return
    df = pd.concat(dfs)

    # convert mutations to int
    df = df.astype({'mutations': 'int32'})
    if 'mutations_low' in df:
        df = df.astype({'mutations_high': 'int32', 'mutations_low': 'int32'})

    # filter mutations by threshold
    df = df[df.mutations > mutations_threshold]

    if df is not None:
        df.to_csv(fname, sep="\t", index=False, float_format='%g')


def read_decomposition(fname):
    """
        Load results of signature decomposition from file (filename)
        Return tuple: numpy array with exposures and a list of signare names
    """
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
