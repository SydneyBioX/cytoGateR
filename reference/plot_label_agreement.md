# Plot agreement between two label columns

Computes per-label agreement (Jaccard overlap) between two label columns
in \`colData(spe)\` and visualises it as a horizontal bar plot with
class sizes.

## Usage

``` r
plot_label_agreement(
  spe,
  label_col1,
  label_col2,
  drop_na = TRUE,
  title = "Label agreement",
  subtitle = "total matches / union"
)
```

## Arguments

- spe:

  A \`SpatialExperiment\`/\`SingleCellExperiment\`.

- label_col1:

  First label column name.

- label_col2:

  Second label column name.

- drop_na:

  Logical; drop rows with NA in either label. Default TRUE.

- title:

  Plot title.

- subtitle:

  Optional subtitle.

## Value

A \`ggplot\` object.

## Examples

``` r
data("cytoGateR_example", package = "cytoGateR")
spe <- cytoGateR_example
plot_label_agreement(spe, "cell_type", "knn_label")
```
