# Introduction to cytoGateR

## 1 Introduction

`cytoGateR` provides marker-aware cell-type annotation for spatial
protein imaging data. It combines interpretable marker gating with
supervised prediction of cells that cannot be confidently gated. The
package operates on `SpatialExperiment` and `SummarizedExperiment`
objects and stores labels, probabilities, and uncertainty measurements
in standard Bioconductor containers.

The package is designed for imaging assays with a relatively small
protein panel, where marker knowledge remains important and annotation
decisions should be inspectable. The main workflow consists of:

1.  defining positive and negative markers for candidate cell types;
2.  estimating cell-type probabilities with tree or soft gating;
3.  selecting high-confidence reference cells;
4.  predicting the remaining cells with a supervised classifier; and
5.  visualising labels, probabilities, and uncertainty.

## 2 Installation

Install the package from Bioconductor:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("cytoGateR")
```

Load the package:

``` r

library(cytoGateR)
```

## 3 Example data

The package includes a compact `SpatialExperiment` derived from the
IMMUcan 2022 Cancer Example dataset. It contains normalized marker
expression, spatial coordinates, tree-gating results, high-confidence
reference labels, and weighted k-nearest-neighbour (kNN) predictions.

``` r

data("cytoGateR_example", package = "cytoGateR")
cytoGateR_example
```

    class: SpatialExperiment
    dim: 40 7204
    metadata(0):
    assays(1): exprs
    rownames(40): MPO H3 ... DNA1 DNA2
    rowData names(0):
    colnames(7204): 12_1 12_2 ... 14_2844 14_2845
    colData names(28): sample_id image_name ... KNN_P_Treg KNN_P_Tumor
    reducedDimNames(0):
    mainExpName: IMMUcan_2022_CancerExample_v1
    altExpNames(0):
    spatialCoords names(2) : cell_x cell_y
    imgData names(1): sample_id

The object contains 40 protein markers:

``` r

head(rownames(cytoGateR_example))
```

    [1] "MPO"    "H3"     "SMA"    "CD16"   "CD38"   "HLA_DR"

``` r

SummarizedExperiment::assayNames(cytoGateR_example)
```

    [1] "exprs"

Sample and image identifiers are stored in `colData()`:

``` r

cd <- SummarizedExperiment::colData(cytoGateR_example)
table(cd$sample_id)
```


    Patient4_006 Patient4_008
            4363         2841 

``` r

table(cd$image_name)
```


    Patient4_006 Patient4_008
            4363         2841 

Spatial coordinates are available through the standard
[`SpatialExperiment::spatialCoords()`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment-methods.html)
accessor:

``` r

head(SpatialExperiment::spatialCoords(cytoGateR_example))
```

           cell_x    cell_y
    12_1  78.0000 1.0857143
    12_2 219.7045 1.0227273
    12_3 302.0800 0.7600000
    12_4 365.7000 0.8500000
    12_5 426.7273 0.4545455
    12_6 473.0909 0.9545455

## 4 Stored annotation outputs

Tree-gating probabilities use the prefix `P_`, while the complete kNN
probability matrix uses the prefix `KNN_P_`.

``` r

tree_probability_columns <- grep("^P_", names(cd), value = TRUE)
knn_probability_columns <- grep("^KNN_P_", names(cd), value = TRUE)

tree_probability_columns
```

     [1] "P_Tumor"       "P_Stroma"      "P_Myeloid"     "P_Neutrophil"
     [5] "P_Bcell"       "P_Plasma_cell" "P_BnTcell"     "P_CD8"
     [9] "P_CD4"         "P_Treg"       

``` r

knn_probability_columns
```

     [1] "KNN_P_Bcell"       "KNN_P_BnTcell"     "KNN_P_CD4"
     [4] "KNN_P_CD8"         "KNN_P_Myeloid"     "KNN_P_Neutrophil"
     [7] "KNN_P_Plasma_cell" "KNN_P_Stroma"      "KNN_P_Treg"
    [10] "KNN_P_Tumor"      

The stored kNN probabilities form a complete probability distribution
for each cell:

``` r

knn_probabilities <- as.matrix(cd[, knn_probability_columns])
summary(rowSums(knn_probabilities))
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
          1       1       1       1       1       1 

Frequently used annotation columns include:

``` r

annotation_columns <- c(
  "cell_type_hard",
  "core_group",
  "knn_label",
  "knn_confidence"
)
annotation_columns[annotation_columns %in% names(cd)]
```

    [1] "cell_type_hard" "core_group"     "knn_label"      "knn_confidence"

The next vignette demonstrates marker-aware tree gating and
visualisation. The third vignette uses high-confidence cells as a
supervised-learning reference.

## 5 Session information

``` r

sessionInfo()
```

    R version 4.6.1 (2026-06-24)
    Platform: x86_64-pc-linux-gnu
    Running under: Ubuntu 24.04.5 LTS

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

    locale:
     [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8
     [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8
     [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C
    [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C

    time zone: UTC
    tzcode source: system (glibc)

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base

    other attached packages:
    [1] cytoGateR_0.99.0

    loaded via a namespace (and not attached):
     [1] SummarizedExperiment_1.42.0 gtable_0.3.6
     [3] rjson_0.2.23                xfun_0.61
     [5] ggplot2_4.0.3               ggrepel_0.9.8
     [7] Biobase_2.72.0              lattice_0.22-9
     [9] vctrs_0.7.3                 tools_4.6.1
    [11] generics_0.1.4              stats4_4.6.1
    [13] proxy_0.4-29                tibble_3.3.1
    [15] pkgconfig_2.0.3             Matrix_1.7-5
    [17] RColorBrewer_1.1-3          S7_0.2.2
    [19] S4Vectors_0.50.3            lifecycle_1.0.5
    [21] compiler_4.6.1              farver_2.1.2
    [23] ggforce_0.5.0               Seqinfo_1.2.0
    [25] graphlayouts_1.2.5          htmltools_0.5.9
    [27] yaml_2.3.12                 pillar_1.11.1
    [29] tidyr_1.3.2                 MASS_7.3-65
    [31] SingleCellExperiment_1.34.0 DelayedArray_0.38.2
    [33] cachem_1.1.0                dbscan_1.2.7
    [35] magick_2.9.1                viridis_0.6.5
    [37] abind_1.4-8                 mclust_6.1.3
    [39] tidyselect_1.2.1            digest_0.6.39
    [41] dplyr_1.2.1                 purrr_1.2.2
    [43] polyclip_1.10-7             fastmap_1.2.0
    [45] grid_4.6.1                  cli_3.6.6
    [47] SparseArray_1.12.3          magrittr_2.0.5
    [49] S4Arrays_1.12.1             dichromat_2.0-1
    [51] ggraph_2.2.2                tidygraph_1.3.1
    [53] withr_3.0.3                 scales_1.4.0
    [55] rmarkdown_2.32              XVector_0.52.0
    [57] matrixStats_1.5.0           igraph_2.3.4
    [59] otel_0.2.0                  gridExtra_2.3.1
    [61] SpatialExperiment_1.22.0    memoise_2.0.1
    [63] evaluate_1.0.5              knitr_1.52
    [65] GenomicRanges_1.64.0        IRanges_2.46.0
    [67] viridisLite_0.4.3           rlang_1.3.0
    [69] Rcpp_1.1.2                  glue_1.8.1
    [71] tweenr_2.0.3                BiocGenerics_0.58.1
    [73] jsonlite_2.0.0              R6_2.6.1
    [75] MatrixGenerics_1.24.0      
