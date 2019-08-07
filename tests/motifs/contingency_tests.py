import pytest
import numpy as np
from mutagene.motifs import *


@pytest.mark.parametrize(
    "contingency_table",
    np.array([
         [
             [20, 10],
             [40, 60]
         ]
     ])
)
def test_OR(contingency_table):
    # purpose: test odds ratio
    contingency_table = Haldane_correction(contingency_table)
    observed = calculate_OR(contingency_table)
    assert observed == 1/3


@pytest.mark.parametrize(
    "contingency_table",
    np.array([
         [
             [20, 10],
             [40, 60]
         ]
     ])
)
def test_RR(contingency_table):
    # purpose: test risk ratio/enrichment
    contingency_table = Haldane_correction(contingency_table)
    observed = calculate_RR(contingency_table)
    assert observed == 5/9


@pytest.mark.parametrize(
    "contingency_table",
    np.array([
         [
             [0, 10],
             [40, 60]
         ]
     ])
)
def test_OR_zero(contingency_table):
    # purpose: test odds ratio
    contingency_table = Haldane_correction(contingency_table)
    observed = calculate_OR(contingency_table)
    assert np.isclose(observed, 14.057851239669422)


@pytest.mark.parametrize(
    "contingency_table",
    np.array([
         [
             [0, 0],
             [40, 60]
         ]
     ])
)
def test_RR_zero(contingency_table):
    # purpose: test risk ratio/enrichment
    contingency_table = Haldane_correction(contingency_table)
    observed = calculate_RR(contingency_table)
    assert np.isclose(observed, 0.8347107438016529)


@pytest.mark.parametrize(
    "contingency_table",
    np.array([
         [
             [20, 10],
             [40, 60]
         ]
     ])
)
def test_bad_stat(contingency_table):
    # purpose: test if user enters incorrect pvalue type
    contingency_table = Haldane_correction(contingency_table)
    observed = get_stats(contingency_table, "XX")
    assert observed == 1.0


@pytest.mark.parametrize(
    "contingency_table",
    np.array([
         [
             [-5, -10],
             [40, 60]
         ]
     ])
)
def test_pval_errors(contingency_table):
    # purpose: test error handling when pvalue error
    contingency_table = Haldane_correction(contingency_table)
    observed = get_stats(contingency_table)
    observed_chi = get_stats(contingency_table, 'chi2')
    assert observed == observed_chi == 1.0


@pytest.mark.parametrize(
    "contingency_table",
    np.array([
         [
             [20, 10],
             [40, 60]
         ]
     ])
)
def test_bad_stat_type(contingency_table):
    # purpose: test risk ratio/enrichment
    contingency_table = Haldane_correction(contingency_table)
    observed = get_stats(contingency_table, "XX")
    assert observed == 1.0


# @pytest.mark.parametrize(
#     "contingency_table",
#     np.array([
#          [
#              [20, 10],
#              [40, 60]
#          ]
#      ])
# )
# def test_pvalue_fisher(contingency_table):
#     # purpose: test Fisher's pvalue
#     observed = get_stats(contingency_table)
#     return observed
#
#
# @pytest.mark.parametrize(
#     "contingency_table",
#     np.array([
#          [
#              [20, 10],
#              [40, 60]
#          ]
#      ])
# )
# def test_pvalue_chi(contingency_table):
#     # purpose: test Chi-Squared pvalue
#     observed = get_stats(contingency_table, "chi2")
#     print(observed)