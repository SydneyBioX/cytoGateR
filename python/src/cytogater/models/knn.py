"""Weighted k-nearest-neighbour annotation."""

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import rankdata

from ..preprocessing import expression_frame


def predict_wknn_multi(
    train_data,
    test_data,
    train_labels,
    k=5,
    method="pearson",
    return_matrix=False,
    chunk_size=10000,
):
    """Port of R ``predict_wknn_multi`` for cells x markers matrices."""
    train = np.asarray(train_data, dtype=float)
    test = np.asarray(test_data, dtype=float)
    labels = np.asarray(train_labels).astype(str)
    if isinstance(train_labels, pd.Series) and isinstance(
        train_labels.dtype, pd.CategoricalDtype
    ):
        classes = np.asarray(train_labels.cat.categories).astype(str)
    elif isinstance(train_labels, pd.Categorical):
        classes = np.asarray(train_labels.categories).astype(str)
    else:
        classes = np.unique(labels)
    if train.ndim != 2 or test.ndim != 2 or train.shape[1] != test.shape[1]:
        raise ValueError(
            "Dimension mismatch: train_data and test_data must be "
            "cells x markers with the same marker columns."
        )
    if len(test) == 0:
        empty = pd.DataFrame(index=getattr(test_data, "index", None), columns=classes)
        result = {"labels": np.array([], dtype=str), "probs": np.array([])}
        if return_matrix:
            result["prob_matrix"] = empty
            result["prob_mat"] = empty
        return result
    if method not in {"pearson", "spearman", "cosine", "euclidean"}:
        raise ValueError(f"Unknown method: {method}")
    chunk_size = len(test) if chunk_size is None or chunk_size <= 0 else int(chunk_size)
    if method == "cosine":
        chunk_size = len(test)
    k = min(int(k), len(train))
    out_labels = np.empty(len(test), dtype=object)
    out_probs = np.zeros(len(test))
    probabilities = np.zeros((len(test), len(classes)))

    train_work = np.apply_along_axis(rankdata, 1, train) if method == "spearman" else train
    for start in range(0, len(test), chunk_size):
        stop = min(start + chunk_size, len(test))
        chunk = test[start:stop]
        chunk_work = (
            np.apply_along_axis(rankdata, 1, chunk)
            if method == "spearman" else chunk
        )
        if method in {"pearson", "spearman"}:
            a = chunk_work - np.nanmean(chunk_work, axis=1, keepdims=True)
            b = train_work - np.nanmean(train_work, axis=1, keepdims=True)
            scores = (a @ b.T) / (
                np.linalg.norm(a, axis=1, keepdims=True)
                * np.linalg.norm(b, axis=1)[None, :]
            )
            is_distance = False
        elif method == "cosine":
            scores = (chunk @ train.T) / (
                np.linalg.norm(chunk, axis=1, keepdims=True)
                * np.linalg.norm(train, axis=1)[None, :]
            )
            is_distance = False
        else:
            scores = cdist(chunk, train, metric="euclidean")
            is_distance = True

        for local_i, row in enumerate(scores):
            global_i = start + local_i
            order = np.argsort(row, kind="stable")
            if not is_distance:
                order = np.argsort(-row, kind="stable")
            neighbors = order[:k]
            weights = (
                1 / (row[neighbors] + 1e-6)
                if is_distance else np.maximum(row[neighbors], 0)
            )
            weights[~np.isfinite(weights)] = 0
            neighbor_labels = labels[neighbors]
            weight_sum = weights.sum()
            if not np.isfinite(weight_sum) or weight_sum <= 0:
                out_labels[global_i] = neighbor_labels[0]
                continue
            label_sums = np.array([
                weights[neighbor_labels == cell_type].sum()
                for cell_type in classes
            ])
            winner = int(np.argmax(label_sums))
            out_labels[global_i] = classes[winner]
            out_probs[global_i] = label_sums[winner] / weight_sum
            probabilities[global_i] = label_sums / weight_sum

    probability_frame = pd.DataFrame(
        probabilities, index=getattr(test_data, "index", None), columns=classes
    )
    result = {
        "labels": out_labels.astype(str),
        "probs": out_probs,
    }
    if return_matrix:
        # R calls this component prob_matrix. Keep prob_mat as a compatibility alias.
        result["prob_matrix"] = probability_frame
        result["prob_mat"] = probability_frame
    return result


def train_custom_knn(
    spe,
    label_col="cutoff_label",
    assay_name="exprs",
    unknown_label="Unknown",
    features="all",
    cv_folds=5,
    repeats=10,
    agreement_thresh=0.8,
    k=5,
    method="pearson",
    seed=None,
    chunk_size=250,
    BPPARAM=None,
):
    """Clean reference labels by repeated CV and return a kNN reference."""
    del BPPARAM
    frame = expression_frame(spe, assay_name)
    features_used = (
        list(frame.columns)
        if isinstance(features, str) and features == "all"
        else [x for x in features if x in frame]
    )
    labels = spe.obs[label_col]
    core = labels.notna() & (labels != unknown_label)
    x, y = frame.loc[core, features_used], labels[core].astype(str)
    class_levels = sorted(y.unique())
    y = pd.Series(
        pd.Categorical(y, categories=class_levels),
        index=y.index,
        name=y.name,
    )
    if len(x) == 0:
        raise ValueError("No labelled core cells are available for kNN training.")
    if cv_folds < 2:
        raise ValueError("cv_folds must be at least 2.")
    agreements = np.zeros(len(x))
    sum_probabilities = np.zeros((len(x), len(class_levels)))
    rng = np.random.default_rng(seed)
    for _ in range(repeats):
        fold_assignment = np.zeros(len(x), dtype=int)
        for cell_type in class_levels:
            class_index = np.flatnonzero(y.to_numpy() == cell_type)
            fold_values = np.resize(np.arange(cv_folds), len(class_index))
            fold_assignment[class_index] = rng.permutation(fold_values)
        for fold in range(cv_folds):
            train = np.flatnonzero(fold_assignment != fold)
            test = np.flatnonzero(fold_assignment == fold)
            if not len(test):
                continue
            pred = predict_wknn_multi(
                x.iloc[train], x.iloc[test], y.iloc[train], k, method,
                return_matrix=True, chunk_size=chunk_size,
            )
            sum_probabilities[test] += pred["prob_matrix"].to_numpy()
            agreements[test] += pred["labels"] == y.iloc[test].astype(str).to_numpy()
    agreements /= repeats
    average_probabilities = sum_probabilities / repeats
    retained = agreements >= agreement_thresh
    cleaned = labels.copy()
    cleaned.loc[x.index[~retained]] = unknown_label
    spe.obs["cleaned_core_label"] = cleaned
    model = {
        "type": "wknn",
        "reference_data": x.loc[retained],
        "reference_labels": y.loc[retained],
        "params": {"k": k, "method": method},
    }
    return {
        "spe": spe,
        "model": model,
        "reference_data": x.loc[retained],
        "reference_labels": y.loc[retained],
        "core_prob_mat": pd.DataFrame(
            average_probabilities[retained],
            index=x.index[retained],
            columns=class_levels,
        ),
        "features_used": features_used,
        "agreement_rates": pd.Series(agreements[retained], index=x.index[retained]),
        "cleaned_core_names": x.index[retained].tolist(),
    }


def predict_unknown_with_knn(
    spe,
    knn_ref,
    assay_name="exprs",
    label_col="cutoff_label",
    out_col="knn_label_filled",
    pred_col="knn_pred",
    unknown_label="Unknown",
    unassigned_label="Unassigned",
    threshold=0.6,
    k=5,
    dist_method="pearson",
    chunk_size=10000,
):
    frame = expression_frame(spe, assay_name)
    current_labels = (
        spe.obs["cleaned_core_label"].copy()
        if "cleaned_core_label" in spe.obs
        else spe.obs[label_col].copy()
    )
    replace = current_labels.isna() | (current_labels == unknown_label)
    if not replace.any():
        return {"spe": spe, "prob_mat": None}
    model = knn_ref.get("model", knn_ref)
    prediction = predict_wknn_multi(
        model["reference_data"],
        frame.loc[replace, knn_ref["features_used"]],
        model["reference_labels"],
        k, dist_method, True, chunk_size,
    )
    probs = prediction["prob_mat"]
    probs.index = frame.index[replace]
    confidence = probs.max(axis=1)
    predicted = probs.idxmax(axis=1).where(confidence >= threshold, unassigned_label)
    spe.obs[pred_col] = pd.Series(index=spe.obs_names, dtype=object)
    spe.obs.loc[replace, pred_col] = predicted
    current_labels.loc[replace] = predicted
    spe.obs[out_col] = current_labels
    spe.obs["knn_confidence"] = pd.Series(confidence, index=confidence.index)
    return {"spe": spe, "prob_mat": probs}


def predict_hierarchical_knn_recursive(
    spe,
    hier_ref,
    hc_tree,
    assay_name="exprs",
    threshold=0.7,
    agreement_threshold=0.8,
    k=5,
    repeats=5,
    dist_methods=("pearson", "cosine"),
    BPPARAM=None,
    out_col="hier_label",
    chunk_size=1000,
    unassigned_label="Unassigned",
):
    """Traverse node references using an ensemble of weighted kNN predictions."""
    del hc_tree, BPPARAM
    frame = expression_frame(spe, assay_name)
    node_by_cluster = {
        value["cluster_id"]: value for key, value in hier_ref.items() if key != "_root"
    }

    def descend(cell, cluster):
        node = node_by_cluster.get(cluster)
        if node is None:
            return unassigned_label
        votes, confidence = [], []
        for method in dist_methods:
            for _ in range(repeats):
                result = predict_wknn_multi(
                    node["train_data"], cell[node["markers"]].to_frame().T,
                    node["train_labels"], k, method, chunk_size=chunk_size,
                )
                votes.append(result["labels"][0])
                confidence.append(result["probs"][0])
        winner = max(set(votes), key=votes.count)
        if np.mean(confidence) < threshold or votes.count(winner) / len(votes) < agreement_threshold:
            return unassigned_label
        members = node["left_members"] if winner == "Left" else node["right_members"]
        if len(members) == 1:
            return members[0]
        candidates = [
            value["cluster_id"] for value in node_by_cluster.values()
            if set(value["left_members"] + value["right_members"]) == set(members)
        ]
        return descend(cell, candidates[0]) if candidates else unassigned_label

    spe.obs[out_col] = [descend(row, hier_ref["_root"]) for _, row in frame.iterrows()]
    return spe
