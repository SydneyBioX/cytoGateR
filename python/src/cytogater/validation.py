"""Minimal validation matching the R package's required inputs."""

import pandas as pd


def validate_lineage_table(lineage_table):
    table = pd.DataFrame(lineage_table)
    missing = {"cell_type", "pos_markers", "neg_markers"} - set(table)
    if missing:
        raise ValueError("lineage_table missing columns: " + ", ".join(sorted(missing)))
    return table


def validate_adata(spe, assay_name="exprs"):
    if assay_name not in (None, "X") and assay_name not in spe.layers:
        raise ValueError(f"Assay '{assay_name}' not found in adata.layers.")
    return spe

