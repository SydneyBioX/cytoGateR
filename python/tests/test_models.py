import numpy as np
import pytest

pytest.importorskip("sklearn")
from cytogater.models.knn import predict_unknown_with_knn, train_custom_knn
from cytogater.models.random_forest import (
    predict_unknown_with_randomforest,
    train_custom_randomforest,
)


def test_knn_train_and_predict(two_lineages):
    adata, _ = two_lineages
    fit = train_custom_knn(
        adata, label_col="core", cv_folds=2, repeats=1,
        agreement_thresh=0, k=3, seed=1,
    )
    result = predict_unknown_with_knn(
        adata, fit, label_col="core", k=3, threshold=0.5,
    )
    assert result["prob_mat"].shape == (30, 2)
    assert adata.obs.loc[adata.obs.core == "Unknown", "knn_label_filled"].notna().all()


def test_random_forest_train_and_predict(two_lineages):
    adata, _ = two_lineages
    fit = train_custom_randomforest(
        adata, label_col="core", num_trees=20, cv_folds=2,
        repeats=1, agreement_thresh=0, seed=1, num_threads=1,
    )
    result = predict_unknown_with_randomforest(
        adata, fit["model"], label_col="core", threshold=0.5,
    )
    assert result["prob_mat"].shape == (30, 2)
    assert np.isfinite(adata.obs.loc[adata.obs.core == "Unknown", "rf_confidence"]).all()


def test_consensus_cleaning_drops_single_cell_classes(two_lineages):
    adata, _ = two_lineages
    adata.obs.loc[adata.obs.index[0], "core"] = "Rare"

    knn = train_custom_knn(
        adata, label_col="core", cv_folds=3, repeats=1,
        agreement_thresh=0.5, k=3, seed=1,
    )
    assert knn["spe"].obs.loc[adata.obs.index[0], "cleaned_core_label"] == "Unknown"

    forest = train_custom_randomforest(
        adata, label_col="core", num_trees=10, cv_folds=3,
        repeats=1, agreement_thresh=0.5, seed=1, num_threads=1,
    )
    assert forest["spe"].obs.loc[adata.obs.index[0], "cleaned_core_label"] == "Unknown"
