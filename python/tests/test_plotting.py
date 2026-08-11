import pytest


matplotlib = pytest.importorskip("matplotlib", reason="requires the [plot] extra")
matplotlib.use("Agg")

from cytogater.plotting import plot_celltype_tree, tree_to_df  # noqa: E402


def _example_tree():
    return {
        "type": "node",
        "depth": 0,
        "cells": list(range(6)),
        "marker": "CD3",
        "cutoff": 0.5,
        "sep_score": 0.8,
        "left": {"type": "leaf", "depth": 1, "cells": [0, 1, 2]},
        "right": {"type": "leaf", "depth": 1, "cells": [3, 4, 5]},
    }


def test_tree_table_matches_r_plot_fields():
    table = tree_to_df(_example_tree())

    assert list(table["id"]) == ["root", "root_L", "root_R"]
    assert table.loc[0, "label"] == "CD3 > 0.50\nsep=0.80"
    assert table.loc[1, "label"] == "leaf\nn=3"
    assert {"parent", "type", "label", "sep_score", "n_cells"}.issubset(table.columns)


def test_celltype_tree_uses_elbow_edges_and_boxed_labels():
    ax = plot_celltype_tree(_example_tree(), title="T cell gating tree")

    assert ax.get_title() == "T cell gating tree"
    assert len(ax.lines) == 2
    assert len(ax.texts) == 3
    assert all(text.get_bbox_patch() is not None for text in ax.texts)
    assert ax.axison is False
