# Confusion matrix between two label columns

Builds a confusion matrix comparing two label columns in
\`colData(spe)\`, with aligned levels so missing categories are
represented.

## Usage

``` r
label_confusion_matrix(spe, label_col1, label_col2, drop_na = TRUE)
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

A table with rows = \`label_col1\` levels and columns = \`label_col2\`
levels.

## Examples

``` r
data("cytoGateR_example", package = "cytoGateR")
label_confusion_matrix(
  cytoGateR_example,
  "cell_type_hard",
  "knn_label"
)
#>              v2
#> v1            Bcell BnTcell  CD4  CD8 Myeloid Neutrophil Plasma_cell Stroma
#>   Bcell         280       7    0    4       5          0           5     16
#>   BnTcell       192    1279  563  375       6          0          12      0
#>   CD4             0       0    5    0       0          0           0      0
#>   CD8             0       0    1  138       0          0           3      2
#>   Myeloid        91       2   81   55     450         67          58    254
#>   Neutrophil      0       0    4    2       3        314           8     53
#>   Plasma_cell     0       0    0    1       0          0          26      4
#>   Stroma          0       0    0    0       0          0           0     80
#>   Treg           10       2  158    1       3          6           8     28
#>   Tumor           0       0    0    0       0          0           0      0
#>   Unassigned      0       0    0    0       0          0           0      0
#>   Uncertain      28      20   79   65       6         19          32    614
#>              v2
#> v1            Treg Tumor Unassigned Uncertain
#>   Bcell          0     0          9         0
#>   BnTcell       87     0         57         0
#>   CD4            0     0          0         0
#>   CD8            0    12          1         0
#>   Myeloid        8     1         42         0
#>   Neutrophil     0   558          2         0
#>   Plasma_cell    0     0          0         0
#>   Stroma         0     0          0         0
#>   Treg         190   215          8         0
#>   Tumor          0    10          0         0
#>   Unassigned     0     0          0         0
#>   Uncertain      7   446         26         0
```
