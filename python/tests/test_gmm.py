import numpy as np

from cytogater.gating.gmm import fit_gmm_2, gmm_equal_posterior_cutoff


def test_gmm_finds_two_components():
    rng = np.random.default_rng(3)
    values = np.r_[rng.normal(0, 0.2, 100), rng.normal(2, 0.2, 100)]
    fit = fit_gmm_2(values)
    assert fit["mu1"] < 0.1 < fit["mu2"]
    assert 0.8 < fit["cutoff"] < 1.2


def test_equal_posterior_symmetric_case():
    assert gmm_equal_posterior_cutoff(0, 2, 1, 1, 0.5, 0.5) == 1

