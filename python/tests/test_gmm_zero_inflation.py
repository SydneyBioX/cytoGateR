"""Zero-inflation behaviour of fit_gmm_2.

These run in milliseconds and need no external data, so they localise the
failure that ``test_r_parity_tree_gating.py`` observes end to end on real
IMC/MIBI panels, where most markers are >75% zero.
"""

import numpy as np
import pytest

from cytogater.core import build_fullcoverage_tree, fit_gmm_2

import pandas as pd


def zero_inflated(zero_fraction, n=5000, positive_mean=1.0, seed=0):
    """A clean two-component signal buried in a spike at zero."""
    rng = np.random.default_rng(seed)
    n_zero = int(round(n * zero_fraction))
    positive = rng.normal(positive_mean, 0.1, n - n_zero)
    return np.concatenate([np.zeros(n_zero), np.clip(positive, 0.05, None)])


def test_separates_a_mildly_zero_inflated_marker():
    """Below the 75% spike the percentile initialisation still straddles zero."""
    fit = fit_gmm_2(zero_inflated(0.50), "equal_posteriors")
    assert fit is not None
    assert fit["sep_score"] > 1.0


@pytest.mark.parametrize("zero_fraction", [0.80, 0.90, 0.95, 0.98, 0.995])
def test_separates_a_heavily_zero_inflated_marker(zero_fraction):
    """The positive mode is unambiguous; only the fitted model used to hide it.

    Regression guard for the equal-variance fit in fit_gmm_2. Two earlier
    versions failed here: seeding at the 25th/75th percentiles put both
    components inside the zero spike, and re-seeding at the 1st/99th percentiles
    fixed that only until the spike passed ~99% and swallowed those quantiles
    too.
    """
    values = zero_inflated(zero_fraction)
    assert values.max() > 0.5, "sanity: the positive mode exists"
    fit = fit_gmm_2(values, "equal_posteriors")
    assert fit is not None
    assert fit["sep_score"] > 1.0


def test_returns_none_when_too_few_cells_carry_the_second_component():
    """5 positive cells in 5,000 is not a two-component signal.

    ``mclust::Mclust(x, G = 2)`` errors out on the same input, so returning None
    here keeps the R and Python versions agreeing about what is unfittable.
    """
    assert fit_gmm_2(zero_inflated(0.999), "equal_posteriors") is None


@pytest.mark.parametrize("zero_fraction", [0.80, 0.95, 0.995])
def test_does_not_collapse_a_component_onto_the_zero_spike(zero_fraction):
    """The real defect: an unequal-variance component shrinking onto the spike.

    Its SD runs to the variance floor, which drives the equal-posterior cutoff
    to ~0 so every non-zero cell is called positive, and leaves a scale-free
    sep_score that still looks plausible.
    """
    values = zero_inflated(zero_fraction)
    fit = fit_gmm_2(values, "equal_posteriors")
    assert fit is not None
    assert min(fit["s1"], fit["s2"]) / max(fit["s1"], fit["s2"]) > 0.05
    assert fit["cutoff"] > 0.01, "cutoff collapsed to the zero spike"
    assert fit["mu1"] < fit["cutoff"] < fit["mu2"]


def test_picks_the_unequal_variance_model_when_it_genuinely_fits_better():
    """Components of clearly different width must not be forced to share one."""
    rng = np.random.default_rng(0)
    values = np.concatenate([rng.normal(0.62, 0.038, 3000), rng.normal(0.78, 0.103, 3000)])
    tied = fit_gmm_2(values, "equal_posteriors", gmm_model_names="E")
    free = fit_gmm_2(values, "equal_posteriors", gmm_model_names="V")
    auto = fit_gmm_2(values, "equal_posteriors")
    assert tied["s1"] == tied["s2"]
    assert free["s1"] != free["s2"]
    assert auto["sep_score"] == pytest.approx(free["sep_score"])


def test_rejects_an_unknown_model_name():
    with pytest.raises(ValueError, match="gmm_model_names"):
        fit_gmm_2(zero_inflated(0.5), "equal_posteriors", gmm_model_names="VVV")


def test_tree_still_splits_on_a_heavily_zero_inflated_marker():
    """A zero sep_score falls below min_score, leaving a bare leaf.

    That leaf is what turned into an all-zero probability column for whole cell
    types on real MIBI panels.
    """
    frame = pd.DataFrame({"CD31": zero_inflated(0.94)})
    tree = build_fullcoverage_tree(
        frame, ["CD31"], [], max_depth=6, min_cells=100, min_score=0.05,
        cutoff_method="equal_posteriors",
    )
    assert tree["type"] == "node"
