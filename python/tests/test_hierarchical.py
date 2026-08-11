import numpy as np

from cytogater.gating.hierarchical import (
    build_hierarchical_reference,
    build_lineage_hierarchy,
)
from cytogater.models import knn as knn_module


def test_hierarchical_prediction_uses_r_style_subsampling_and_chunks(
    monkeypatch, two_lineages
):
    adata, _ = two_lineages
    adata.obs["cleaned_core_label"] = adata.obs["core"]
    tree = build_lineage_hierarchy(adata)
    reference = build_hierarchical_reference(
        adata, tree, unknown_label="Unknown", top_n=None
    )

    calls = []
    original = knn_module.predict_wknn_multi

    def recording_predict(train_data, test_data, train_labels, *args, **kwargs):
        calls.append((len(train_data), len(test_data)))
        return original(train_data, test_data, train_labels, *args, **kwargs)

    monkeypatch.setattr(knn_module, "predict_wknn_multi", recording_predict)
    np.random.seed(11)
    knn_module.predict_hierarchical_knn_recursive(
        adata,
        reference,
        tree,
        threshold=0,
        agreement_threshold=0,
        k=3,
        repeats=2,
        dist_methods=("euclidean",),
        chunk_size=25,
        out_col="hierarchical_label",
    )

    root = reference["Node_1"]
    expected_training_size = int(np.floor(0.8 * len(root["train_data"])))
    assert calls
    assert all(training_size == expected_training_size for training_size, _ in calls)
    assert max(test_size for _, test_size in calls) > 1
    assert max(test_size for _, test_size in calls) <= 25
    assert set(adata.obs["hierarchical_label"]) <= set(tree["labels"])
    assert adata.obs["hierarchical_label"].notna().all()
