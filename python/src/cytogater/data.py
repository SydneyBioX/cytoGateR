"""Small data conversion helpers."""

import pandas as pd


def lineage_table(cell_type, pos_markers, neg_markers=None):
    """Create the list-column lineage table expected by cytoGateR."""
    if neg_markers is None:
        neg_markers = [[] for _ in cell_type]
    return pd.DataFrame({
        "cell_type": cell_type,
        "pos_markers": pos_markers,
        "neg_markers": neg_markers,
    })


def probabilities_from_obs(spe, prefix="P_"):
    """Read probability columns stored in ``adata.obs``."""
    columns = [x for x in spe.obs.columns if x.startswith(prefix)]
    out = spe.obs[columns].copy()
    out.columns = [x[len(prefix):].replace("_", " ") for x in columns]
    return out

