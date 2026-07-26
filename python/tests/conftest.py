import numpy as np
import pandas as pd
import pytest


class MiniAnnData:
    def __init__(self, x, markers):
        self.X = x
        self.layers = {"exprs": x}
        self.var_names = pd.Index(markers)
        self.obs_names = pd.Index([f"cell_{i}" for i in range(len(x))])
        self.obs = pd.DataFrame(index=self.obs_names)
        self.obsm = {}
        self.uns = {}


@pytest.fixture
def two_lineages():
    rng = np.random.default_rng(1)
    low = rng.normal(0, 0.15, (60, 2))
    high = rng.normal(3, 0.15, (60, 2))
    x = np.vstack([
        np.column_stack([high[:, 0], low[:, 0]]),
        np.column_stack([low[:, 1], high[:, 1]]),
    ])
    adata = MiniAnnData(x, ["CD3", "CD20"])
    adata.obs["core"] = ["Tcell"] * 45 + ["Unknown"] * 15 + ["Bcell"] * 45 + ["Unknown"] * 15
    table = pd.DataFrame({
        "cell_type": ["Tcell", "Bcell"],
        "pos_markers": [["CD3"], ["CD20"]],
        "neg_markers": [["CD20"], ["CD3"]],
    })
    return adata, table

