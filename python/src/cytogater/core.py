"""Compact Python port of the cytoGateR soft- and tree-gating core."""

from __future__ import annotations

import re
from math import log, pi

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import expit


def _expression(spe, assay_name):
    """Return cells x markers expression and marker names from AnnData."""
    x = spe.X if assay_name in (None, "X") else spe.layers[assay_name]
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x, dtype=float), np.asarray(spe.var_names, dtype=str)


def _lineage(lineage_table, markers):
    table = pd.DataFrame(lineage_table).copy()
    available = set(markers)
    table["cell_type"] = table["cell_type"].astype(str).str.strip()
    table["pos_markers"] = table["pos_markers"].map(
        lambda values: [str(x) for x in values if str(x) in available]
    )
    table["neg_markers"] = table["neg_markers"].map(
        lambda values: [str(x) for x in values if str(x) in available]
    )
    return table[table["pos_markers"].map(bool)].reset_index(drop=True)


def gmm_equal_posterior_cutoff(mu1, mu2, s1, s2, p1, p2):
    """Cutoff between two Gaussian components with equal posterior density."""
    values = np.asarray([mu1, mu2, s1, s2, p1, p2], dtype=float)
    midpoint = float(np.mean([mu1, mu2]))
    if not np.isfinite(values).all() or s1 <= 0 or s2 <= 0:
        return midpoint

    def difference(x):
        a = p1 * np.exp(-0.5 * ((x - mu1) / s1) ** 2) / s1
        b = p2 * np.exp(-0.5 * ((x - mu2) / s2) ** 2) / s2
        return a - b

    try:
        return float(brentq(difference, min(mu1, mu2), max(mu1, mu2)))
    except ValueError:
        return midpoint


def fit_gmm_2(x, cutoff_method="mean", gmm_model_names=None):
    """Fit the R package's two-component, one-dimensional GMM."""
    del gmm_model_names  # mclust-specific; retained for API alignment.
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 50 or len(np.unique(x)) < 3:
        return None

    # Small deterministic 1-D EM avoids exposing extra Python-only parameters.
    mu = np.quantile(x, [0.25, 0.75]).astype(float)
    sd = np.repeat(max(float(np.std(x)), 1e-6), 2)
    weight = np.array([0.5, 0.5])
    previous = -np.inf
    for _ in range(200):
        density = np.column_stack(
            [
                weight[k] / (sd[k] * np.sqrt(2 * pi))
                * np.exp(-0.5 * ((x - mu[k]) / sd[k]) ** 2)
                for k in range(2)
            ]
        )
        total = density.sum(axis=1, keepdims=True)
        total[total == 0] = np.finfo(float).tiny
        responsibility = density / total
        nk = responsibility.sum(axis=0)
        mu = (responsibility * x[:, None]).sum(axis=0) / nk
        variance = (
            responsibility * (x[:, None] - mu) ** 2
        ).sum(axis=0) / nk
        sd = np.sqrt(np.maximum(variance, 1e-12))
        weight = nk / len(x)
        likelihood = float(np.log(total).sum())
        if abs(likelihood - previous) < 1e-6:
            break
        previous = likelihood

    order = np.argsort(mu)
    mu1, mu2 = mu[order]
    s1, s2 = sd[order]
    p1, p2 = weight[order]
    cutoff = (
        gmm_equal_posterior_cutoff(mu1, mu2, s1, s2, p1, p2)
        if cutoff_method == "equal_posteriors"
        else float((mu1 + mu2) / 2)
    )
    return {
        "mu1": float(mu1), "mu2": float(mu2),
        "s1": float(s1), "s2": float(s2),
        "p1": float(p1), "p2": float(p2),
        "cutoff": cutoff,
        "sep_score": float(abs(mu2 - mu1) / np.sqrt(s1**2 + s2**2)),
    }


def marker_separability(x):
    fit = fit_gmm_2(x)
    return np.nan if fit is None else fit["sep_score"]


def score_marker_logistic(x, cutoff, scale):
    return expit((np.asarray(x, dtype=float) - cutoff) / scale)


def score_marker_rank(x):
    values = pd.Series(np.asarray(x, dtype=float))
    ranks = values.rank(method="average", na_option="keep").to_numpy()
    return (ranks - 1) / (np.isfinite(ranks).sum() - 1 + 1e-9)


def fit_marker_stats(spe, markers, assay_name="exprs", min_n=50):
    x, marker_names = _expression(spe, assay_name)
    positions = {name: i for i, name in enumerate(marker_names)}
    result = {}
    for marker in markers:
        if marker not in positions:
            continue
        values = x[:, positions[marker]]
        values = values[np.isfinite(values)]
        mad = np.median(np.abs(values - np.median(values))) * 1.4826
        mad = mad if np.isfinite(mad) and mad != 0 else 1e-3
        fallback = {
            "valid": False,
            "cutoff": float(np.quantile(values, 0.98)) if len(values) else 0.0,
            "scale": float(mad),
            "weight": 0.1,
            "mu_high": float(np.max(values)) if len(values) else np.nan,
        }
        fit = fit_gmm_2(values) if len(values) >= min_n else None
        if fit is None:
            result[marker] = fallback
        else:
            result[marker] = {
                "valid": True,
                "cutoff": (fit["mu1"] + fit["mu2"]) / 2,
                "scale": (fit["s1"] + fit["s2"]) / 2 or mad,
                "weight": fit["sep_score"] if fit["sep_score"] > 0 else 0.1,
                "mu_high": fit["mu2"],
            }
    return result


def assign_soft_labels(prob_mat, unknown_thresh=0.4):
    values = np.asarray(prob_mat)
    names = np.asarray(
        prob_mat.columns if hasattr(prob_mat, "columns")
        else [str(i) for i in range(values.shape[1])]
    )
    labels = names[np.argmax(values, axis=1)].astype(object)
    labels[np.max(values, axis=1) < unknown_thresh] = "Unknown"
    return labels


def calculate_soft_scores(
    spe, marker_stats, lineage_table, assay_name="exprs", neg_strength=0.8
):
    x, marker_names = _expression(spe, assay_name)
    positions = {name: i for i, name in enumerate(marker_names)}
    table = _lineage(lineage_table, marker_names)
    scores = np.zeros((x.shape[0], len(table)))
    for j, row in table.iterrows():
        weighted, weight_sum = np.zeros(x.shape[0]), 0.0
        for marker in row.pos_markers:
            stat = marker_stats[marker]
            weight = stat["weight"] if np.isfinite(stat["weight"]) else 0.1
            probability = score_marker_logistic(
                x[:, positions[marker]], stat["cutoff"], stat["scale"]
            )
            probability[~np.isfinite(probability)] = 0
            weighted += probability * weight
            weight_sum += weight
        base = weighted / weight_sum if weight_sum > 0 else weighted
        penalty = np.ones(x.shape[0])
        for marker in row.neg_markers:
            stat = marker_stats[marker]
            bad = score_marker_logistic(
                x[:, positions[marker]], stat["cutoff"], stat["scale"]
            )
            bad[~np.isfinite(bad)] = 0
            penalty *= 1 - bad * neg_strength
        scores[:, j] = base * penalty
    return pd.DataFrame(scores, index=spe.obs_names, columns=table.cell_type)


def build_fullcoverage_tree(
    expr_mat,
    markers_pos,
    markers_neg=(),
    cell_idx=None,
    depth=0,
    max_depth=4,
    min_cells=200,
    min_score=0.5,
    cutoff_method="mean",
    gmm_model_names=None,
):
    """Build the same recursive marker-priority tree as cytoGateR R."""
    frame = expr_mat if isinstance(expr_mat, pd.DataFrame) else pd.DataFrame(expr_mat)
    markers_pos = [m for m in markers_pos if m in frame.columns]
    markers_neg = [m for m in markers_neg if m in frame.columns]
    cell_idx = np.arange(len(frame)) if cell_idx is None else np.asarray(cell_idx)
    leaf = {"type": "leaf", "depth": depth, "cells": cell_idx}
    if len(cell_idx) < min_cells or depth >= max_depth or not markers_pos:
        return leaf

    subset = frame.iloc[cell_idx]
    pos_score = subset[markers_pos].mean(axis=1).to_numpy()
    neg_score = subset[markers_neg].mean(axis=1).to_numpy() if markers_neg else 0
    target = (pos_score - neg_score > np.nanmedian(pos_score - neg_score)).astype(int)
    candidates = []
    for marker in markers_pos:
        values = subset[marker].to_numpy(float)
        fit = fit_gmm_2(values, cutoff_method, gmm_model_names)
        if fit is None:
            continue
        scale = np.sqrt(fit["s1"] ** 2 + fit["s2"] ** 2)
        predicted = (values > fit["cutoff"]).astype(int)
        accuracy = max(np.mean(predicted == target), np.mean(1 - predicted == target))
        candidates.append((fit["sep_score"] + 2 * (accuracy - 0.5), marker, fit, scale))
    if not candidates:
        return leaf
    _, marker, fit, scale = max(candidates, key=lambda item: item[0])
    if not np.isfinite(fit["sep_score"]) or fit["sep_score"] < min_score:
        return leaf
    values = frame.iloc[cell_idx][marker].to_numpy(float)
    left, right = cell_idx[values <= fit["cutoff"]], cell_idx[values > fit["cutoff"]]
    if not len(left) or not len(right):
        return leaf
    remaining = [m for m in markers_pos if m != marker]
    return {
        "type": "node", "depth": depth, "marker": marker,
        "cutoff": fit["cutoff"], "scale": float(scale),
        "sep_score": fit["sep_score"], "cells": cell_idx,
        "left": build_fullcoverage_tree(
            frame, remaining, markers_neg, left, depth + 1, max_depth,
            min_cells, min_score, cutoff_method, gmm_model_names
        ),
        "right": build_fullcoverage_tree(
            frame, remaining, markers_neg, right, depth + 1, max_depth,
            min_cells, min_score, cutoff_method, gmm_model_names
        ),
    }


def collect_path_scores(tree, expr_mat, cell_i):
    frame = expr_mat if isinstance(expr_mat, pd.DataFrame) else pd.DataFrame(expr_mat)
    scores, node = [], tree
    while node is not None and node["type"] == "node":
        value = float(frame.iloc[cell_i][node["marker"]])
        scores.append(
            float(score_marker_logistic(value, node["cutoff"], node["scale"]))
            if np.isfinite(value) else np.nan
        )
        node = node["right"] if np.isfinite(value) and value > node["cutoff"] else node["left"]
    return np.asarray(scores)


def neg_penalty(
    expr_mat, neg_markers, cell_i, marker_stats=None, neg_strength=0.8
):
    frame = expr_mat if isinstance(expr_mat, pd.DataFrame) else pd.DataFrame(expr_mat)
    marker_stats = marker_stats or {}
    penalty = 1.0
    for marker in set(neg_markers).intersection(marker_stats):
        value = float(frame.iloc[cell_i][marker])
        stat = marker_stats[marker]
        if np.isfinite(value):
            penalty *= 1 - float(
                score_marker_logistic(value, stat["cutoff"], stat["scale"])
            ) * neg_strength
    return float(np.clip(penalty, 0, 1))


def tree_prob(
    tree,
    expr_mat,
    cell_i,
    combine="mean",
    neg_markers=(),
    marker_stats=None,
    neg_strength=0.8,
    lambda_=2,
):
    scores = collect_path_scores(tree, expr_mat, cell_i)
    scores = scores[np.isfinite(scores)]
    if not len(scores):
        return np.nan
    base = np.mean(scores) if combine == "mean" else np.prod(scores)
    return float(
        base * len(scores) / (len(scores) + lambda_)
        * neg_penalty(expr_mat, neg_markers, cell_i, marker_stats, neg_strength)
    )


def run_soft_gating(
    spe, lineage_table, assay_name="exprs", unknown_thresh=0.4, store=True
):
    x, markers = _expression(spe, assay_name)
    del x
    table = _lineage(lineage_table, markers)
    all_markers = list(dict.fromkeys(sum(table.pos_markers.tolist(), []) + sum(table.neg_markers.tolist(), [])))
    stats = fit_marker_stats(spe, all_markers, assay_name)
    probabilities = calculate_soft_scores(spe, stats, table, assay_name)
    probabilities[~np.isfinite(probabilities)] = 0
    labels = assign_soft_labels(probabilities, unknown_thresh)
    if store:
        for cell_type in probabilities:
            spe.obs["P_" + re.sub(r"\s+", "_", cell_type)] = probabilities[cell_type].to_numpy()
        spe.obs["soft_tree_label"] = labels
        spe.obsm["cytogater_soft_probabilities"] = probabilities.to_numpy()
        spe.uns["cytogater_soft_cell_types"] = probabilities.columns.tolist()
    return {
        "spe": spe, "marker_stats": stats, "prob_mat": probabilities,
        "labels": labels, "lineage_table": table,
    }


def run_tree_gating(
    spe,
    lineage_table,
    assay_name="exprs",
    max_depth=4,
    min_cells=200,
    min_score=0.5,
    uncert_thresh=0.25,
    neg_strength=0.4,
    cutoff_method="mean",
    gmm_model_names=None,
    workers=1,
):
    """Run tree gating; parameters and defaults align with the R function."""
    del workers  # Retained for API alignment; compact implementation is serial.
    x, markers = _expression(spe, assay_name)
    frame = pd.DataFrame(x, index=spe.obs_names, columns=markers)
    table = _lineage(lineage_table, markers)
    all_markers = list(dict.fromkeys(sum(table.pos_markers.tolist(), []) + sum(table.neg_markers.tolist(), [])))
    stats = fit_marker_stats(spe, all_markers, assay_name)
    trees = {
        row.cell_type: build_fullcoverage_tree(
            frame, row.pos_markers, row.neg_markers, max_depth=max_depth,
            min_cells=min_cells, min_score=min_score,
            cutoff_method=cutoff_method, gmm_model_names=gmm_model_names,
        )
        for row in table.itertuples()
    }
    probabilities = pd.DataFrame(0.0, index=spe.obs_names, columns=list(trees))
    for row in table.itertuples():
        probabilities[row.cell_type] = [
            tree_prob(
                trees[row.cell_type], frame, i, "mean", row.neg_markers,
                stats, neg_strength,
            )
            for i in range(len(frame))
        ]
    probabilities[~np.isfinite(probabilities)] = 0
    lineage = [x for x in probabilities.columns if x != "Proliferating"]
    if lineage:
        best_label = probabilities[lineage].idxmax(axis=1).to_numpy(object)
        best_probability = probabilities[lineage].max(axis=1).to_numpy()
    else:
        best_label = np.repeat("Uncertain", len(frame)).astype(object)
        best_probability = np.zeros(len(frame))
    hard_label = best_label.copy()
    hard_label[best_probability < uncert_thresh] = "Uncertain"
    proliferating = (
        probabilities["Proliferating"].to_numpy()
        if "Proliferating" in probabilities else np.zeros(len(frame))
    )
    state = np.where(
        (proliferating >= 0.5) & (hard_label != "Uncertain"),
        np.char.add("Prolif ", hard_label.astype(str)),
        hard_label,
    )
    spe.obs["prob_best"] = best_probability
    spe.obs["cell_type_hard"] = pd.Categorical(hard_label)
    spe.obs["cell_type_state"] = pd.Categorical(state)
    spe.obs["prob_prolif"] = proliferating
    for cell_type in probabilities:
        spe.obs["P_" + re.sub(r"\s+", "_", cell_type)] = probabilities[cell_type].to_numpy()
    spe.obsm["cytogater_tree_probabilities"] = probabilities.to_numpy()
    spe.uns["cytogater_tree_cell_types"] = probabilities.columns.tolist()
    return {
        "spe": spe, "trees": trees, "prob_mat": probabilities,
        "hard_label": hard_label, "lineage_table": table,
    }

