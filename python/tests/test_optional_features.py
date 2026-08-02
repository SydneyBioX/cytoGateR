import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib", reason="requires the [plot] extra")
matplotlib.use("Agg")
pytest.importorskip("torch", reason="requires the [torch] extra")
ad = pytest.importorskip("anndata")

from cytogater.models.neural_network import (  # noqa: E402
    predict_unknown_with_dl,
    train_custom_dl,
)
from cytogater.plotting import plot_labelled_cells  # noqa: E402
from cytogater.uncertainty import calculate_uncertainty  # noqa: E402


def test_anndata_torch_plotting_and_uncertainty():
    rng = np.random.default_rng(7)
    x = np.r_[
        rng.normal([3, 0], 0.2, (30, 2)),
        rng.normal([0, 3], 0.2, (30, 2)),
    ]
    adata = ad.AnnData(x)
    adata.var_names = ["CD3", "CD20"]
    adata.obs_names = [f"cell_{i}" for i in range(60)]
    adata.layers["exprs"] = x.copy()
    adata.obs["label"] = (
        ["T"] * 25 + ["Unknown"] * 5 + ["B"] * 25 + ["Unknown"] * 5
    )
    adata.obs["sample_id"] = "sample_1"
    adata.obsm["spatial"] = rng.normal(size=(60, 2))

    fit = train_custom_dl(
        adata, label_col="label", epochs=1, batch_size=16,
        hidden_dims=(8,), seed=1,
    )
    prediction = predict_unknown_with_dl(
        adata, fit["model"], label_col="label"
    )
    assert prediction["prob_mat"].shape == (10, 2)

    probabilities = pd.DataFrame(
        np.tile([0.8, 0.2], (60, 1)),
        index=adata.obs_names,
        columns=["T", "B"],
    )
    uncertainty = calculate_uncertainty(
        probabilities, adata, k_spatial=3
    )
    assert uncertainty.shape == (60, 5)
    assert plot_labelled_cells(adata, "soft_tree_label_filled") is not None

