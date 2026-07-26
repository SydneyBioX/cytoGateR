import numpy as np
import pandas as pd

import cytogater


class MiniAnnData:
    """Only the AnnData fields used by the gating core."""

    def __init__(self, x, markers):
        self.X = x
        self.layers = {"exprs": x}
        self.var_names = pd.Index(markers)
        self.obs_names = pd.Index([f"cell_{i}" for i in range(len(x))])
        self.obs = pd.DataFrame(index=self.obs_names)
        self.obsm = {}
        self.uns = {}


def example_data():
    rng = np.random.default_rng(1)
    low = rng.normal(0, 0.15, (60, 2))
    high = rng.normal(3, 0.15, (60, 2))
    x = np.vstack([
        np.column_stack([high[:, 0], low[:, 0]]),
        np.column_stack([low[:, 1], high[:, 1]]),
    ])
    table = pd.DataFrame({
        "cell_type": ["Tcell", "Bcell"],
        "pos_markers": [["CD3"], ["CD20"]],
        "neg_markers": [["CD20"], ["CD3"]],
    })
    return MiniAnnData(x, ["CD3", "CD20"]), table


def test_scoring_matches_r_definitions():
    np.testing.assert_allclose(
        cytogater.score_marker_logistic([0, 1, 2], cutoff=1, scale=0.5),
        [0.11920292, 0.5, 0.88079708],
    )
    ranks = cytogater.score_marker_rank([3, 1, 2, np.nan])
    np.testing.assert_allclose(ranks[:3], [1, 0, 0.5], atol=1e-9)
    assert np.isnan(ranks[3])


def test_gmm_and_equal_posterior_cutoff():
    rng = np.random.default_rng(4)
    values = np.r_[rng.normal(0, 0.2, 100), rng.normal(2, 0.2, 100)]
    fit = cytogater.fit_gmm_2(values)
    assert fit["mu1"] < 0.1
    assert fit["mu2"] > 1.9
    assert 0.8 < fit["cutoff"] < 1.2
    assert cytogater.gmm_equal_posterior_cutoff(0, 2, 1, 1, 0.5, 0.5) == 1


def test_soft_gating_pipeline_stores_r_style_outputs():
    adata, table = example_data()
    result = cytogater.run_soft_gating(adata, table)
    assert (result["labels"][:60] == "Tcell").mean() > 0.95
    assert (result["labels"][60:] == "Bcell").mean() > 0.95
    assert list(result["prob_mat"].columns) == ["Tcell", "Bcell"]
    assert "P_Tcell" in adata.obs
    assert "soft_tree_label" in adata.obs


def test_tree_gating_pipeline_stores_r_style_outputs():
    adata, table = example_data()
    result = cytogater.run_tree_gating(
        adata,
        table,
        max_depth=1,
        min_cells=20,
        min_score=0.1,
        uncert_thresh=0.1,
    )
    assert list(result["trees"]) == ["Tcell", "Bcell"]
    assert result["prob_mat"].shape == (120, 2)
    assert (result["hard_label"][:60] == "Tcell").mean() > 0.95
    assert (result["hard_label"][60:] == "Bcell").mean() > 0.95
    assert "cell_type_hard" in adata.obs
    assert "cytogater_tree_probabilities" in adata.obsm

