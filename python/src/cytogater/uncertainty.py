"""Probability and spatial uncertainty."""

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


def calculate_uncertainty(prob_mat, spe, sample_col="sample_id", k_spatial=15):
    probabilities = pd.DataFrame(prob_mat).reindex(spe.obs_names)
    values = probabilities.to_numpy(float)
    labels = probabilities.columns[
        np.nanargmax(np.nan_to_num(values, nan=-np.inf), axis=1)
    ].to_numpy()
    log_values = np.zeros_like(values)
    positive = values > 0
    log_values[positive] = np.log(values[positive])
    entropy = -np.nansum(values * log_values, axis=1) / np.log(values.shape[1])
    gini = (1 - np.nansum(values**2, axis=1)) / (1 - 1 / values.shape[1])
    ordered = np.sort(values, axis=1)[:, ::-1]
    margin = 1 - (ordered[:, 0] - ordered[:, 1])
    discordance = np.full(len(values), np.nan)
    coordinates = np.asarray(spe.obsm["spatial"])
    for sample in pd.unique(spe.obs[sample_col]):
        index = np.flatnonzero(spe.obs[sample_col].to_numpy() == sample)
        if len(index) <= k_spatial:
            continue
        neighbors = NearestNeighbors(n_neighbors=k_spatial + 1).fit(
            coordinates[index]
        ).kneighbors(return_distance=False)[:, 1:]
        discordance[index] = np.mean(labels[index][neighbors] != labels[index, None], axis=1)
    return pd.DataFrame({
        "cell_id": spe.obs_names, "entropy": entropy,
        "gini_impurity": gini, "margin_uncertainty": margin,
        "spatial_discordance": discordance,
    }, index=spe.obs_names)


def calculate_spatial_prior_labels(
    spe,
    prob_mat,
    k_spatial=50,
    lambda_=0.2,
    protect_threshold=0.85,
    out_col="knn_spatial_label",
):
    probabilities = pd.DataFrame(prob_mat).reindex(spe.obs_names)
    values = probabilities.to_numpy(float)
    classes = probabilities.columns.to_numpy()
    winner = np.nanargmax(np.nan_to_num(values, nan=-np.inf), axis=1)
    confidence = np.nanmax(values, axis=1)
    k = min(k_spatial, len(values) - 1)
    neighbors = NearestNeighbors(n_neighbors=k + 1).fit(
        spe.obsm["spatial"]
    ).kneighbors(return_distance=False)[:, 1:]
    result = classes[winner].astype(object)
    for i in range(len(values)):
        if confidence[i] >= protect_threshold:
            continue
        prior = np.bincount(winner[neighbors[i]], minlength=len(classes)) / k
        posterior = values[i] * prior**lambda_
        result[i] = classes[np.nanargmax(np.nan_to_num(posterior))]
    spe.obs[out_col] = result
    return spe
