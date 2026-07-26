from cytogater.gating.soft import run_soft_gating


def test_soft_gating(two_lineages):
    adata, table = two_lineages
    result = run_soft_gating(adata, table)
    assert (result["labels"][:60] == "Tcell").mean() > 0.95
    assert (result["labels"][60:] == "Bcell").mean() > 0.95
    assert "soft_tree_label" in adata.obs

