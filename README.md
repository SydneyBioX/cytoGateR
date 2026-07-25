# cytoGateR

`cytoGateR` is an R package for marker-aware cell-type annotation in spatial
protein imaging data, such as imaging mass cytometry (IMC) and CODEX data.
It combines interpretable marker gating with supervised prediction to annotate
cells that cannot be confidently assigned by gating alone.

The package works with standard Bioconductor `SpatialExperiment` and
`SummarizedExperiment` objects. Annotation labels, class probabilities, and
uncertainty measurements are stored in the object metadata so that results
remain compatible with other Bioconductor workflows.

## Development status

`cytoGateR` is under active development and is being prepared for submission to
Bioconductor. The `main` branch is the development branch and may include
experimental material that is not part of the final Bioconductor software
package.

A separate package-only branch or repository will be used for Bioconductor
submission and maintenance.

## Main features

- Marker-aware soft and tree gating
- Positive- and negative-marker definitions for each cell type
- High-confidence reference-cell selection
- Weighted k-nearest-neighbour prediction
- Random-forest classification
- Optional Torch-based deep learning
- Recursive hierarchical kNN classification
- Probability, confidence, and spatial uncertainty measurements
- Visualization of marker distributions, labels, probabilities, trees, and
  model agreement

## Installation

### Development version

Install the development version from GitHub:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  "SydneyBioX/cytoGateR",
  dependencies = TRUE
)
```

The optional deep-learning functions additionally require the R `torch`
package and its runtime:

```r
install.packages("torch")
torch::install_torch()
```

### Bioconductor version

After the package is accepted into Bioconductor, the released version will be
installed with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("cytoGateR")
```

## Quick start

Load the package and its included spatial proteomics example:

```r
library(cytoGateR)

data("cytoGateR_example", package = "cytoGateR")
cytoGateR_example
```

Select one complete image:

```r
image_ids <- SummarizedExperiment::colData(cytoGateR_example)$image_name
spe <- cytoGateR_example[, image_ids == image_ids[1L]]
```

Define a small marker lineage table:

```r
lineage_table <- data.frame(
  cell_type = c("Tcell", "Bcell"),
  pos_markers = I(list("CD3e", "CD20")),
  neg_markers = I(list("CD20", "CD3e"))
)
```

Run tree gating:

```r
tree_result <- run_tree_gating(
  spe,
  lineage_table,
  assay_name = "exprs",
  max_depth = 1,
  min_cells = 100,
  min_score = 0.1,
  workers = 1
)

table(tree_result$hard_label)
```

Select high-confidence reference cells:

```r
core_result <- apply_cutoff_labels(
  tree_result,
  cutoff_fn = function(x) {
    stats::quantile(x, 0.95, na.rm = TRUE)
  },
  label_col = "core_group",
  unknown_label = "Unknown"
)
```

The packaged demo also contains full tree-gating and weighted-kNN results:

```r
cd <- SummarizedExperiment::colData(cytoGateR_example)

table(cd$core_group)
table(cd$knn_label)

knn_probability_columns <- grep(
  "^KNN_P_",
  names(cd),
  value = TRUE
)

summary(rowSums(as.matrix(cd[, knn_probability_columns])))
```

## Vignettes

The package includes four evaluated vignettes:

1. **Introduction to cytoGateR** — package concepts and example-data structure
2. **Tree Gating and Visualisation** — marker-aware gating and spatial plots
3. **Machine-Learning Prediction** — weighted kNN training, prediction, and
   uncertainty
4. **Hierarchical kNN Classification** — lineage hierarchy construction and
   recursive prediction

After installation, list the available vignettes with:

```r
browseVignettes("cytoGateR")
```

## Recreating the example data

The packaged example is derived from the IMMUcan 2022 Cancer Example dataset
provided by `imcdatasets`. Its complete creation workflow is documented in:

```text
inst/scripts/create_cytoGateR_example.R
```

From the R package source directory, regenerate it with:

```bash
Rscript inst/scripts/create_cytoGateR_example.R
```

The script performs sample selection, metadata reduction, tree gating,
high-confidence reference selection, weighted-kNN prediction, probability
validation, and compressed `.rda` creation.

## Development

Run the unit tests with:

```r
devtools::test()
```

Regenerate documentation with:

```r
roxygen2::roxygenise()
```

Build the package and its Quarto vignettes with:

```bash
R CMD build cytoGateR
```

Issues and feature requests can be reported through the
[GitHub issue tracker](https://github.com/SydneyBioX/cytoGateR/issues).

## License

`cytoGateR` is distributed under the GPL-2 license.
