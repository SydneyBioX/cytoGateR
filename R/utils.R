#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @keywords internal
.safe_mad <- function(x, na.rm = TRUE) {
  v <- stats::mad(x, na.rm = na.rm)
  if (!is.finite(v) || v == 0) 1e-6 else v
}

#' @keywords internal
.assert_spe <- function(spe) {
  ok <- inherits(spe, c("SpatialExperiment", "SingleCellExperiment"))
  if (!ok) stop("spe must be a SpatialExperiment or SingleCellExperiment.")
}

#' @keywords internal
.assert_lineage_table <- function(lineage_table) {
  req <- c("cell_type", "pos_markers", "neg_markers")
  miss <- setdiff(req, names(lineage_table))
  if (length(miss) > 0) stop("lineage_table missing columns: ", paste(miss, collapse = ", "))

  # minimal type checks
  if (!is.list(lineage_table$pos_markers) || !is.list(lineage_table$neg_markers)) {
    stop("lineage_table$pos_markers and neg_markers must be list-columns (each element a character vector).")
  }
}

#' @keywords internal
.clean_lineage_table <- function(lineage_table, spe) {

  lineage_table <- dplyr::as_tibble(lineage_table)

  lineage_table |>
    dplyr::mutate(
      cell_type   = trimws(as.character(.data$cell_type)),
      pos_markers = lapply(.data$pos_markers, function(v) intersect(unlist(v), rownames(spe))),
      neg_markers = lapply(.data$neg_markers, function(v) intersect(unlist(v), rownames(spe)))
    ) |>
    dplyr::filter(lengths(.data$pos_markers) > 0)
}
