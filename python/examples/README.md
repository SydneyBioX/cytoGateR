# Examples

## Basic notebook

`basic_pipeline.ipynb` shows the shortest tree-gating workflow.

## Complete R-data pipeline

`full_pipeline.ipynb` loads the packaged R `cytoGateR_example.rda`
`SpatialExperiment`, converts one complete image to `AnnData`, and runs:

1. soft gating;
2. tree gating;
3. high-confidence reference selection;
4. weighted kNN consensus cleaning and prediction;
5. random-forest consensus cleaning and prediction;
6. hierarchical kNN;
7. uncertainty and protected spatial priors;
8. metrics, plots, probability exports, and AnnData output.

Start Jupyter from the repository and open the notebook:

```bash
conda activate cytogater-python
jupyter lab examples/full_pipeline.ipynb
```

The final cell saves the converted, annotated object as
`cytogater_r_demo_python.h5ad`.
