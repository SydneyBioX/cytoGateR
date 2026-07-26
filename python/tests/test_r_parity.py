from pathlib import Path

import pandas as pd
import numpy as np

from cytogater import score_marker_logistic
from cytogater.models.knn import predict_wknn_multi


def test_logistic_scores_match_r_fixture():
    fixture = pd.read_csv(
        Path(__file__).parent / "fixtures" / "r_scoring_reference.csv"
    )
    actual = score_marker_logistic(
        fixture.value, fixture.cutoff, fixture.scale
    )
    np.testing.assert_allclose(actual, fixture.r_logistic, rtol=1e-12)


def test_weighted_pearson_knn_matches_r_reference():
    train = pd.DataFrame(
        [[1, 2, 3], [3, 2, 1], [1, 1, 0]],
        index=["t1", "t2", "t3"],
    )
    test = pd.DataFrame(
        [[1, 2, 2.8], [2.8, 2, 1]],
        index=["q1", "q2"],
    )
    labels = pd.Series(
        pd.Categorical(["A", "B", "C"], categories=["A", "B", "C"]),
        index=train.index,
    )
    result = predict_wknn_multi(
        train, test, labels, k=2, method="pearson", return_matrix=True
    )
    np.testing.assert_array_equal(result["labels"], ["A", "B"])
    np.testing.assert_allclose(result["probs"], [1, 0.5268425], atol=1e-7)
    np.testing.assert_allclose(
        result["prob_matrix"],
        [[1, 0, 0], [0, 0.5268425, 0.4731575]],
        atol=1e-7,
    )
