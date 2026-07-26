"""Result types shared by applications using cytoGateR."""

from typing import Any, TypedDict

import numpy as np
import pandas as pd


class GatingResult(TypedDict):
    spe: Any
    prob_mat: pd.DataFrame
    lineage_table: pd.DataFrame
    hard_label: np.ndarray
    trees: dict[str, dict]


class PredictionResult(TypedDict):
    spe: Any
    prob_mat: pd.DataFrame | None

