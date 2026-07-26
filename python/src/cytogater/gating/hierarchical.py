"""Hierarchical reference construction."""

import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform

from ..preprocessing import expression_frame


def build_lineage_hierarchy(
    spe, label_col="cleaned_core_label", assay_name="exprs"
):
    """Build a complete-linkage hierarchy from cell-type pseudobulks."""
    frame = expression_frame(spe, assay_name)
    labels = spe.obs[label_col]
    keep = labels.notna() & (labels.astype(str) != "Unknown")
    groups = list(dict.fromkeys(labels[keep].astype(str)))
    averages = np.column_stack([
        frame.loc[keep & (labels.astype(str) == group)].mean(axis=0)
        for group in groups
    ])
    correlation = np.corrcoef(averages, rowvar=False)
    distance = np.nan_to_num(1 - correlation, nan=1.0)
    distance = np.clip((distance + distance.T) / 2, 0, None)
    np.fill_diagonal(distance, 0)
    return {"linkage": linkage(squareform(distance), method="complete"), "labels": groups}


def build_hierarchical_reference(
    spe,
    hc_tree,
    marker_stats=None,
    label_col="cleaned_core_label",
    top_n=5,
    assay_name="exprs",
    unknown_label="Unassigned",
):
    """Create the binary training reference for every hierarchy node."""
    frame = expression_frame(spe, assay_name)
    labels = spe.obs[label_col].astype(object)
    keep = labels.notna() & (labels != unknown_label)
    root, nodes = to_tree(hc_tree["linkage"], rd=True)
    leaf_names = hc_tree["labels"]

    def leaves(node):
        if node.is_leaf():
            return [leaf_names[node.id]]
        return leaves(node.left) + leaves(node.right)

    result = {}
    internal = [node for node in nodes if not node.is_leaf()]
    for number, node in enumerate(internal, 1):
        left_members, right_members = leaves(node.left), leaves(node.right)
        selected = keep & labels.isin(left_members + right_members)
        node_labels = np.where(labels[selected].isin(left_members), "Left", "Right")
        markers = list(marker_stats) if marker_stats is not None else list(frame.columns)
        if top_n is not None and marker_stats is not None:
            left_mean = frame.loc[selected].loc[node_labels == "Left", markers].mean()
            right_mean = frame.loc[selected].loc[node_labels == "Right", markers].mean()
            score = (left_mean - right_mean).abs() * np.array(
                [marker_stats[m].get("weight", 0.1) for m in markers]
            )
            markers = score.sort_values(ascending=False).index[:top_n].tolist()
        result[f"Node_{number}"] = {
            "train_data": frame.loc[selected, markers],
            "train_labels": node_labels,
            "markers": markers,
            "left_members": left_members,
            "right_members": right_members,
            "cluster_id": node.id,
        }
    result["_root"] = root.id
    return result
