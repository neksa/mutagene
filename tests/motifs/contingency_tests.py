import numpy as np

from mutagene.motifs import *


def test_OR():
    """purpose: test odds ratio"""
    contingency_table = make_contingency_table(np.array([[20, 10], [40, 60]]))
    contingency_table = Haldane_correction(contingency_table)
    observed = calculate_OR(contingency_table)
    assert observed == 3.0


def test_RR():
    """purpose: test risk ratio/enrichment"""
    contingency_table = make_contingency_table(np.array([[20, 10], [40, 60]]))
    contingency_table = Haldane_correction(contingency_table)
    observed = calculate_RR(contingency_table)
    assert np.isclose(observed, 5 / 3)


def test_OR_zero():
    """purpose: test odds ratio with zero"""
    contingency_table = make_contingency_table(np.array([[0, 10], [40, 60]]))
    contingency_table = Haldane_correction(contingency_table)
    observed = calculate_OR(contingency_table)
    assert np.isclose(observed, 0.07113462669018224)


def test_RR_zero():
    """purpose: test risk ratio/enrichment with zeros"""
    contingency_table = make_contingency_table(np.array([[0, 0], [40, 60]]))
    contingency_table = Haldane_correction(contingency_table)
    observed = calculate_RR(contingency_table)
    assert np.isclose(observed, 1.2469135802469136)


def test_bad_stat():
    """purpose: test if user enters incorrect pvalue type"""
    contingency_table = make_contingency_table(np.array([[20, 10], [40, 60]]))
    contingency_table = Haldane_correction(contingency_table)
    observed = get_stats(contingency_table, "XX")
    assert observed == 1.0


def test_pval_errors():
    """purpose: test error handling when pvalue error"""
    contingency_table = make_contingency_table(np.array([[-5, -10], [40, 60]]))
    contingency_table = Haldane_correction(contingency_table)
    observed = get_stats(contingency_table)
    observed_chi = get_stats(contingency_table, "chi2")
    assert observed == observed_chi == 1.0


def test_bad_stat_type():
    """purpose: test risk ratio/enrichment"""
    contingency_table = make_contingency_table(np.array([[20, 10], [40, 60]]))
    contingency_table = Haldane_correction(contingency_table)
    observed = get_stats(contingency_table, "XX")
    assert observed == 1.0
