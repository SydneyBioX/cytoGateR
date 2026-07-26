"""Expression access and marker preprocessing."""

import numpy as np
import pandas as pd


def expression_frame(spe, assay_name="exprs"):
    x = spe.X if assay_name in (None, "X") else spe.layers[assay_name]
    if hasattr(x, "toarray"):
        x = x.toarray()
    return pd.DataFrame(
        np.asarray(x, dtype=float), index=spe.obs_names, columns=spe.var_names
    )


def clean_lineage_table(lineage_table, spe):
    table = pd.DataFrame(lineage_table).copy()
    markers = set(map(str, spe.var_names))
    table["cell_type"] = table.cell_type.astype(str).str.strip()
    for column in ("pos_markers", "neg_markers"):
        table[column] = table[column].map(
            lambda values: [str(x) for x in values if str(x) in markers]
        )
    return table[table.pos_markers.map(bool)].reset_index(drop=True)

