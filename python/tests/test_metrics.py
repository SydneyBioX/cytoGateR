import numpy as np
import pandas as pd

from cytogater.metrics import (
    apply_cutoff_labels,
    assign_confident_labels,
    class_metrics_from_fit,
    probability_label_matrix,
)


def test_class_metrics_match_expected_counts():
    table = class_metrics_from_fit(
        {}, ["B", "B", "T", "T"], ["B", "T", "T", "T"]
    ).set_index("class")
    assert table.loc["B", "tp"] == 1
    assert table.loc["B", "fn"] == 1
    assert table.loc["T", "fp"] == 1


def test_probability_labels_and_confidence():
    probabilities = pd.DataFrame(
        [[0.9, 0.1], [0.5, 0.5]], columns=["B", "T"]
    )
    labels = assign_confident_labels(probabilities, [True, False])
    np.testing.assert_array_equal(labels, ["B", "Unknown"])
    matrix = probability_label_matrix(probabilities, lambda x: 0.8)
    assert matrix.iloc[0].tolist() == [True, False]


def test_cutoff_labels_require_exactly_one_positive(two_lineages):
    adata, _ = two_lineages
    probabilities = pd.DataFrame(
        [[0.9, 0.1], [0.9, 0.9], [0.1, 0.1]],
        index=adata.obs_names[:3],
        columns=["B", "T"],
    )
    # Use a three-cell AnnData-like view for alignment with the result matrix.
    adata.obs = adata.obs.iloc[:3].copy()
    adata.obs_names = adata.obs.index
    result = apply_cutoff_labels(
        {"spe": adata, "prob_mat": probabilities},
        cutoff_fn=lambda column: 0.8,
    )
    assert result["spe"].obs["cutoff_label"].tolist() == [
        "B", "Unknown", "Unknown"
    ]
