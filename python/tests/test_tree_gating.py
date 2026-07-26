from cytogater.gating.tree import run_tree_gating


def test_tree_gating(two_lineages):
    adata, table = two_lineages
    result = run_tree_gating(
        adata, table, max_depth=1, min_cells=20,
        min_score=0.1, uncert_thresh=0.1,
    )
    assert result["prob_mat"].shape == (120, 2)
    assert (result["hard_label"][:60] == "Tcell").mean() > 0.95
    assert (result["hard_label"][60:] == "Bcell").mean() > 0.95

