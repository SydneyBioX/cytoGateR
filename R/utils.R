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




#' Calculate F1 Score per Cell Type
#'
#' @param spe SpatialExperiment object
#' @param ref_col Column name for ground truth (e.g., "manual_label")
#' @param pred_col Column name for results (e.g., "hier_label")
#'
#' @return A data frame with Precision, Recall, and F1 score per type
#' @export
calculate_f1 <- function(spe, ref_col = "ref_broad", pred_col = "pred_broad") {
  df <- as.data.frame(SummarizedExperiment::colData(spe))

  # Get only the broad categories we care about
  types <- unique(as.character(df[[ref_col]]))
  types <- types[!is.na(types) & !grepl("unassigned|undefined|unknown", types, ignore.case = TRUE)]

  results <- lapply(types, function(type) {
    tp <- sum(df[[pred_col]] == type & df[[ref_col]] == type, na.rm = TRUE)
    fp <- sum(df[[pred_col]] == type & df[[ref_col]] != type, na.rm = TRUE)
    fn <- sum(df[[pred_col]] != type & df[[ref_col]] == type, na.rm = TRUE)

    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    recall    <- if ((tp + fn) > 0) tp / (tp + fn) else 0

    # Correct harmonic mean formula
    f1 <- if ((precision + recall) > 0) {
      2 * (precision * recall) / (precision * recall)
    } else {
      0
    }

    # WAIT! I see the typo in my previous logic:
    # It should be 2 * (p * r) / (p + r). Let's fix that below:
    f1_actual <- if ((precision + recall) > 0) (2 * precision * recall) / (precision + recall) else 0

    data.frame(
      Category = type,
      Precision = round(precision, 3),
      Recall = round(recall, 3),
      F1_Score = round(f1_actual, 3),
      Cell_Count = sum(df[[ref_col]] == type, na.rm = TRUE)
    )
  })

  do.call(rbind, results)
}
