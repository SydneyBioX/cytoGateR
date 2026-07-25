#' Minimal spatial proteomics example
#'
#' A small `SpatialExperiment` containing selected samples for demonstrating
#' cytoGateR visualization and annotation functions.
#'
#' @format A `SpatialExperiment` with 40 markers and the selected cells.
#' The object contains:
#' \describe{
#'   \item{assay}{An `"exprs"` expression assay.}
#'   \item{colData}{Sample identifiers, image names, reference labels,
#'     tree-gating labels and probabilities, high-confidence core labels, and
#'     weighted k-nearest-neighbour predictions, confidence values, and full
#'     `KNN_P_*` class-probability columns.}
#'   \item{spatialCoords}{Two-dimensional cell coordinates.}
#' }
#'
#' @source A subset of the IMMUcan 2022 Cancer Example dataset from
#' the `imcdatasets` package.
#'
#' @usage data(cytoGateR_example)
#' @keywords datasets
"cytoGateR_example"
