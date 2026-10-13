# Run tree-based soft gating on a SpatialExperiment

Build a marker-based, hierarchical soft-gating model for each lineage in
a marker table, compute per-cell membership probabilities, and derive
hard lineage/state labels.

## Usage

``` r
run_tree_gating(
  spe,
  lineage_table,
  assay_name = "exprs",
  max_depth = 4,
  min_cells = 200,
  min_score = 0.5,
  uncert_thresh = 0.25,
  neg_strength = 0.4,
  cutoff_method = c("mean", "equal_posteriors"),
  gmm_model_names = NULL,
  workers = 1
)
```

## Arguments

- spe:

  A `SpatialExperiment` or `SingleCellExperiment` containing per-cell
  marker expression.

- lineage_table:

  A `data.frame`/tibble with at least the columns `cell_type`,
  `pos_markers`, and `neg_markers`. Marker columns must be list-columns
  of character vectors. Markers not present in `spe` are dropped during
  cleaning.

- assay_name:

  Character scalar naming the assay in `spe` used as the expression
  matrix (default `"exprs"`).

- max_depth:

  Integer; maximum recursion depth of each gating tree.

- min_cells:

  Integer; minimum number of cells required to allow a node split.

- min_score:

  Numeric; minimum separability score required to accept a split.

- uncert_thresh:

  Numeric in \[0,1\]; cells whose best (non-proliferating) lineage
  probability is below this threshold are labeled `"Uncertain"`.

- neg_strength:

  Numeric in \[0,1\]; strength of the negative-marker penalty applied
  within
  [`tree_prob()`](https://sydneybiox.github.io/cytoGateR/reference/tree_prob.md)
  (default `0.4`). Larger values increase the penalty from
  negative-marker expression.

- cutoff_method:

  Cutoff selection method passed to
  [`fit_gmm_2()`](https://sydneybiox.github.io/cytoGateR/reference/fit_gmm_2.md).
  One of `"mean"` or `"equal_posteriors"`.

- gmm_model_names:

  Optional character vector of model names to pass to
  [`mclust::Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.html)
  for 1D GMM fitting (e.g., `"V"`, `"E"`). Use `NULL` (default) to allow
  selection by BIC.

- workers:

  Optional integer number of workers for parallel execution. Defaults to
  1.

## Value

A named list with components:

- spe:

  The input `spe` with added `colData` columns: `P_*` (per-lineage
  probabilities), `cell_type_hard`, `cell_type_state`, `prob_best`, and
  `prob_prolif`.

- trees:

  A named list of fitted gating trees, one per `cell_type`.

- prob_mat:

  Numeric matrix of per-cell probabilities (cells \\\times\\ cell
  types).

- hard_label:

  Character vector of hard lineage labels after applying
  `uncert_thresh`.

- lineage_table:

  The cleaned marker table actually used to build the trees.

## Details

For each `cell_type` in `lineage_table`, a gating tree is fit from
positive and (optionally) negative marker sets. Per-cell probabilities
are then computed by evaluating each fitted tree via
[`tree_prob()`](https://sydneybiox.github.io/cytoGateR/reference/tree_prob.md),
with an optional negative-marker penalty controlled by `neg_strength`.

Hard labels are assigned as the lineage with the highest probability
(excluding `"Proliferating"`), unless the best probability is below
`uncert_thresh`, in which case the cell is labeled `"Uncertain"`. If a
`"Proliferating"` lineage is present, cells with proliferating
probability \\\ge 0.5\\ and not `"Uncertain"` receive a `"Prolif "`
prefix in `cell_type_state`.

Marker statistics used for logistic scoring/penalization are precomputed
once across all markers appearing in `lineage_table` via
[`fit_marker_stats()`](https://sydneybiox.github.io/cytoGateR/reference/fit_marker_stats.md).
Trees are built with
[`build_fullcoverage_tree()`](https://sydneybiox.github.io/cytoGateR/reference/build_fullcoverage_tree.md)
and evaluated with
[`tree_prob()`](https://sydneybiox.github.io/cytoGateR/reference/tree_prob.md).

## Examples

``` r
### This is a minimal example, for the full example, please see vignettes
data("cytoGateR_example", package = "cytoGateR")
image_ids <- SummarizedExperiment::colData(cytoGateR_example)$image_name
spe <- cytoGateR_example[, image_ids == image_ids[1L]]
lineage_table <- data.frame(
  cell_type = c("Tcell", "Bcell"),
  pos_markers = I(list("CD3e", "CD20")),
  neg_markers = I(list("CD20", "CD3e"))
)
tree_result <- run_tree_gating(
  spe,
  lineage_table,
  max_depth = 1,
  min_cells = 100,
  min_score = 0.1,
  workers = 1
)
#> Building 1th tree - Tcell Cell Type.
#> finished tree 1 - Tcell | elapsed: 0.053 s | peak RAM: 0.02 GB
#> Building 2th tree - Bcell Cell Type.
#> finished tree 2 - Bcell | elapsed: 0.033 s | peak RAM: 0.02 GB
head(tree_result$hard_label)
#> [1] "Uncertain" "Uncertain" "Uncertain" "Uncertain" "Uncertain" "Tcell"    
```
