#' Build Hierarchical Reference with Flexible Marker Selection
#'
#' @param spe SpatialExperiment containing core cells.
#' @param hc_tree The hclust object from build_lineage_hierarchy.
#' @param marker_stats Output from fit_marker_stats(). Required if top_n is not NULL.
#' @param label_col Column with core labels (e.g., "cleaned_core_label").
#' @param top_n Number of top DE markers to use per node. Set to NULL to use all markers.
#' @param assay_name Assay to use (default "exprs").
#' @export
build_hierarchical_reference <- function(spe,
                                         hc_tree,
                                         marker_stats = NULL,
                                         label_col = "cleaned_core_label",
                                         top_n = 5,
                                         assay_name = "exprs") {

  .assert_spe(spe)

  # 1. Extract core data
  cd <- as.data.frame(SummarizedExperiment::colData(spe))
  core_idx <- which(!is.na(cd[[label_col]]) & cd[[label_col]] != "Unknown")
  core_spe <- spe[, core_idx]

  n_nodes <- nrow(hc_tree$merge)
  node_list <- list()

  # Helper to find all leaf indices under a node
  get_leaves <- function(side, hc) {
    if (side < 0) return(hc$labels[-side])
    row <- hc$merge[side, ]
    return(c(get_leaves(row[1], hc), get_leaves(row[2], hc)))
  }

  # 2. Iterate through internal nodes
  for (i in seq_len(n_nodes)) {
    left_labels <- get_leaves(hc_tree$merge[i, 1], hc_tree)
    right_labels <- get_leaves(hc_tree$merge[i, 2], hc_tree)

    # Subset cells for this specific split
    node_cells <- core_spe[, core_spe[[label_col]] %in% c(left_labels, right_labels)]
    node_labels <- ifelse(node_cells[[label_col]] %in% left_labels, "Left", "Right")
    expr_sub <- SummarizedExperiment::assay(node_cells, assay_name)

    # 3. Marker Selection Logic
    all_markers <- if(!is.null(marker_stats)) names(marker_stats) else rownames(spe)

    if (!is.null(top_n) && !is.null(marker_stats)) {
      # NODE-SPECIFIC DE LOGIC
      de_scores <- vapply(all_markers, function(m) {
        vals <- as.numeric(expr_sub[m, ])
        m_left <- mean(vals[node_labels == "Left"], na.rm = TRUE)
        m_right <- mean(vals[node_labels == "Right"], na.rm = TRUE)
        diff <- abs(m_left - m_right)
        # Weight by global quality (weight from marker_stats)
        return(diff * (marker_stats[[m]]$weight %||% 0.1))
      }, numeric(1))

      # Select only the most discriminatory markers for this branch
      relevant_markers <- all_markers[order(de_scores, decreasing = TRUE)][1:min(top_n, length(all_markers))]
    } else {
      # ALL MARKERS LOGIC (Consensus style)
      relevant_markers <- all_markers
    }

    # 4. Store training data
    node_list[[paste0("Node_", i)]] <- list(
      train_data = as.data.frame(t(expr_sub[relevant_markers, , drop = FALSE])),
      train_labels = factor(node_labels),
      markers = relevant_markers,
      left_members = left_labels,
      right_members = right_labels
    )
  }

  return(node_list)
}





#' Build Lineage Hierarchy from Core Cells
#'
#' @param spe A SpatialExperiment object with gated core labels.
#' @param label_col The column containing gated labels (e.g., "cleaned_core_label").
#' @param assay_name Assay to use for pseudobulk calculation.
#'
#' @return An hclust object representing the cell type hierarchy.
#' @export
build_lineage_hierarchy <- function(spe,
                                    label_col = "cleaned_core_label",
                                    assay_name = "exprs") {

  .assert_spe(spe)
  cd <- as.data.frame(SummarizedExperiment::colData(spe))

  # 1. Filter to core cells only (exclude Unknown)
  core_cells <- spe[, !is.na(cd[[label_col]]) & cd[[label_col]] != "Unknown"]

  # 2. Calculate Pseudobulk (Mean expression per cell type)
  feat_mat <- SummarizedExperiment::assay(core_cells, assay_name)
  labels <- SummarizedExperiment::colData(core_cells)[[label_col]]

  # Compute column means for each cell type
  avg_expr <- sapply(unique(labels), function(ct) {
    rowMeans(feat_mat[, labels == ct, drop = FALSE])
  })

  # 3. Hierarchical Clustering (The 'Tree')
  # We use correlation distance as it's more robust for cell types
  dist_mat <- as.dist(1 - cor(avg_expr))
  hc <- stats::hclust(dist_mat, method = "complete")

  return(hc)
}
