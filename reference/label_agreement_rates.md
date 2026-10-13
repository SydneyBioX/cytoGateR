# Agreement rates between two label columns

Computes per-label agreement between two label columns in
\`colData(spe)\`. For each label (union of both columns), agreement is
defined as the number of cells where both labels match that value
divided by the number of cells where either column uses that value
(Jaccard-style overlap).

## Usage

``` r
label_agreement_rates(spe, label_col1, label_col2, drop_na = TRUE)
```

## Arguments

- spe:

  A \`SpatialExperiment\`/\`SingleCellExperiment\`.

- label_col1:

  First label column name.

- label_col2:

  Second label column name.

- drop_na:

  Logical; if TRUE, rows with NA in either label are dropped. Default
  TRUE.

## Value

A data frame with columns \`label\`, \`match_n\`, \`union_n\`, and
\`agreement\`.

## Examples

``` r
data("cytoGateR_example", package = "cytoGateR")
label_agreement_rates(
  cytoGateR_example,
  "cell_type_hard",
  "knn_label"
)
#>          label match_n union_n   agreement
#> 1        Bcell     280     647 0.432766615
#> 2      BnTcell    1279    2602 0.491544965
#> 3          CD4       5     891 0.005611672
#> 4          CD8     138     660 0.209090909
#> 5      Myeloid     450    1132 0.397526502
#> 6   Neutrophil     314    1036 0.303088803
#> 7  Plasma_cell      26     157 0.165605096
#> 8       Stroma      80    1051 0.076117983
#> 9         Treg     190     731 0.259917921
#> 10       Tumor      10    1242 0.008051530
#> 11  Unassigned       0     145 0.000000000
#> 12   Uncertain       0    1342 0.000000000
```
