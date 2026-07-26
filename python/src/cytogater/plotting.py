"""Matplotlib visualisations corresponding to the R plotting API."""

import numpy as np
import pandas as pd

from .preprocessing import expression_frame


def _plt():
    import matplotlib.pyplot as plt
    return plt


def plot_probability_map(
    spe, probability, title=None, spatial_key="spatial", ax=None, **kwargs
):
    ax = ax or _plt().subplots()[1]
    coordinates = np.asarray(spe.obsm[spatial_key])
    values = spe.obs[probability] if isinstance(probability, str) else probability
    artist = ax.scatter(coordinates[:, 0], coordinates[:, 1], c=values, **kwargs)
    ax.set(title=title or str(probability), aspect="equal")
    _plt().colorbar(artist, ax=ax)
    return ax


def plot_labelled_cells(
    spe, label_col, spatial_key="spatial", ax=None, title=None, **kwargs
):
    ax = ax or _plt().subplots()[1]
    coordinates = np.asarray(spe.obsm[spatial_key])
    labels = pd.Categorical(spe.obs[label_col])
    ax.scatter(coordinates[:, 0], coordinates[:, 1], c=labels.codes, **kwargs)
    ax.set(title=title or label_col, aspect="equal")
    return ax


def plot_marker_density(spe, marker, assay_name="exprs", label_col=None, ax=None):
    ax = ax or _plt().subplots()[1]
    frame = expression_frame(spe, assay_name)
    if label_col is None:
        ax.hist(frame[marker].dropna(), bins=40, density=True, alpha=0.7)
    else:
        for label, index in spe.obs.groupby(label_col, observed=True).groups.items():
            ax.hist(frame.loc[index, marker].dropna(), bins=40, density=True,
                    histtype="step", label=str(label))
        ax.legend()
    ax.set(xlabel=marker, ylabel="Density")
    return ax


def plot_confusion_matrix(conf_mat, title="Confusion matrix", ax=None):
    ax = ax or _plt().subplots()[1]
    matrix = np.asarray(conf_mat)
    image = ax.imshow(matrix, cmap="Blues")
    labels = getattr(conf_mat, "columns", range(matrix.shape[1]))
    rows = getattr(conf_mat, "index", range(matrix.shape[0]))
    ax.set(xticks=range(len(labels)), xticklabels=labels,
           yticks=range(len(rows)), yticklabels=rows, title=title)
    _plt().colorbar(image, ax=ax)
    return ax


def plot_class_metrics(class_metrics, title="Per-class precision/recall/F1", subtitle=None, ax=None):
    del subtitle
    ax = ax or _plt().subplots()[1]
    table = pd.DataFrame(class_metrics).set_index("class")
    table[["precision", "recall", "f1"]].plot.bar(ax=ax)
    ax.set(title=title, ylim=(0, 1), ylabel="Score")
    return ax


def plot_label_counts(spe, label_col, ax=None):
    ax = ax or _plt().subplots()[1]
    spe.obs[label_col].value_counts().plot.bar(ax=ax)
    ax.set(ylabel="Cells", title=label_col)
    return ax


def plot_label_agreement(agreement, ax=None):
    ax = ax or _plt().subplots()[1]
    table = pd.DataFrame(agreement)
    ax.bar(table["label"], table["agreement"])
    ax.set(ylim=(0, 1), ylabel="Agreement")
    return ax


def plot_label_confusion_matrix(spe, label_col1, label_col2, ax=None):
    from .metrics import label_confusion_matrix
    return plot_confusion_matrix(
        label_confusion_matrix(spe, label_col1, label_col2), ax=ax
    )


def plot_score_hist(res, ax=None):
    ax = ax or _plt().subplots()[1]
    pd.DataFrame(res["prob_mat"]).plot.hist(alpha=0.4, bins=30, ax=ax)
    return ax


def plot_label_cardinality(label_mat, ax=None):
    ax = ax or _plt().subplots()[1]
    pd.DataFrame(label_mat).sum(axis=1).value_counts().sort_index().plot.bar(ax=ax)
    ax.set(xlabel="Labels per cell", ylabel="Cells")
    return ax


def plot_ct_marker_intensity(spe, markers=None, label_col="cell_type_hard", assay_name="exprs", ax=None):
    ax = ax or _plt().subplots()[1]
    frame = expression_frame(spe, assay_name)
    markers = list(frame.columns) if markers is None else markers
    means = frame[markers].groupby(spe.obs[label_col], observed=True).mean()
    image = ax.imshow(means, aspect="auto", cmap="viridis")
    ax.set(xticks=range(len(markers)), xticklabels=markers,
           yticks=range(len(means)), yticklabels=means.index)
    _plt().colorbar(image, ax=ax)
    return ax


plot_pseudobulk_heatmap = plot_ct_marker_intensity
plot_label_dotplot = plot_ct_marker_intensity


def tree_to_df(node, parent=None, id="root"):
    rows = [{"id": id, "parent": parent, "type": node["type"],
             "marker": node.get("marker"), "cutoff": node.get("cutoff")}]
    if node["type"] == "node":
        return pd.concat(
            [
                pd.DataFrame(rows),
                tree_to_df(node["left"], id, id + "L"),
                tree_to_df(node["right"], id, id + "R"),
            ],
            ignore_index=True,
        )
    return pd.DataFrame(rows)


def print_celltype_tree(node, indent=""):
    text = indent + (
        f"{node['marker']} > {node['cutoff']:.3g}"
        if node["type"] == "node" else "leaf"
    )
    print(text)
    if node["type"] == "node":
        print_celltype_tree(node["left"], indent + "  ")
        print_celltype_tree(node["right"], indent + "  ")


def plot_celltype_tree(tree, title="", ax=None):
    ax = ax or _plt().subplots()[1]
    table = tree_to_df(tree)
    depth = {row.id: len(row.id) - 4 for row in table.itertuples()}
    for row in table.itertuples():
        if row.parent is not None:
            ax.plot([depth[row.parent], depth[row.id]], [table.index[table.id == row.parent][0], row.Index], "k-")
        ax.text(depth[row.id], row.Index, row.marker or "leaf")
    ax.set(title=title)
    ax.axis("off")
    return ax


def plot_marker_priority_tree(lineage_table, marker_stats, cell_type_name, ax=None):
    row = pd.DataFrame(lineage_table).query("cell_type == @cell_type_name").iloc[0]
    values = [marker_stats[x]["weight"] for x in row.pos_markers]
    ax = ax or _plt().subplots()[1]
    ax.bar(row.pos_markers, values)
    ax.set(title=cell_type_name, ylabel="Separability weight")
    return ax


def plot_rand_cell_probs(spe=None, prob_mat=None, n=20, seed=None, ax=None):
    del spe
    frame = pd.DataFrame(prob_mat)
    sample = frame.sample(min(n, len(frame)), random_state=seed)
    ax = ax or _plt().subplots()[1]
    sample.plot.bar(stacked=True, ax=ax)
    return ax
