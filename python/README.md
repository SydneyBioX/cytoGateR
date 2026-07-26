# cytoGateR Python

This is the Python implementation of `cytoGateR`. It uses an `AnnData` object
and keeps the R function names, parameter names, and defaults where Python
syntax permits.

Implemented modules include:

- marker GMMs, soft gating, and tree gating;
- hierarchical reference construction and recursive hierarchical kNN;
- weighted kNN consensus cleaning and prediction;
- random-forest consensus cleaning and prediction;
- optional PyTorch training and prediction;
- confidence labels, cutoffs, F1, confusion, and agreement metrics;
- probability and spatial uncertainty;
- Matplotlib spatial, marker, label, metric, and tree plots.

## Install

```bash
pip install -e cytoGateR/python
```

For an isolated development and test environment:

```bash
cd cytoGateR/python
conda env create -f environment.yml
conda activate cytogater-python
pytest
```

Plotting and neural-network support are optional:

```bash
pip install -e 'cytoGateR/python[plot]'
pip install -e 'cytoGateR/python[torch]'
pip install -e 'cytoGateR/python[notebook]'
```

## Data layout

Python `AnnData` stores cells in rows and markers in columns. Put the expression
matrix used by the R `assay_name = "exprs"` default in:

```python
adata.layers["exprs"]
```

The lineage table is a pandas data frame with `cell_type`, `pos_markers`, and
`neg_markers` columns. The marker columns contain Python lists.

## Example

```python
import pandas as pd
import cytogater

lineage_table = pd.DataFrame({
    "cell_type": ["Tcell", "Bcell"],
    "pos_markers": [["CD3e"], ["CD20"]],
    "neg_markers": [["CD20"], ["CD3e"]],
})

result = cytogater.run_tree_gating(
    adata,
    lineage_table,
    assay_name="exprs",
    max_depth=1,
    min_cells=100,
    min_score=0.1,
    workers=1,
)

result["hard_label"]
result["prob_mat"]
```

Results are returned with the same main keys as R. Labels and `P_*` columns are
also stored in `adata.obs`; probability matrices are stored in `adata.obsm`.

The package is intentionally flat at its public surface:

```python
cytogater.run_soft_gating(...)
cytogater.run_tree_gating(...)
cytogater.train_custom_knn(...)
cytogater.train_custom_randomforest(...)
cytogater.calculate_uncertainty(...)
```

The implementations are also grouped under `cytogater.gating`,
`cytogater.models`, `cytogater.metrics`, `cytogater.uncertainty`, and
`cytogater.plotting`.

## R compatibility note

The Python implementation uses a deterministic two-component Gaussian-mixture
fit. R uses `mclust`, so fitted cutoffs can differ slightly even though the
pipeline, parameters, defaults, scoring, and output structure are aligned.
The Python `lambda_` argument corresponds to the R `lambda` argument because
`lambda` is reserved syntax in Python. `num_trees` corresponds to R's
`num.trees`; `train_custom_randomforest()` also accepts `**{"num.trees": 200}`.

Run tests with:

```bash
cd cytoGateR/python
pytest
```
