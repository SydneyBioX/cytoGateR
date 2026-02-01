#' Plot a marker priority "tree" for a given cell type (weights order)
#'
#' This is a simple visualization: markers ordered by weight and connected
#' as a chain ROOT -> marker1 -> marker2 -> ...
#'
#' @param lineage_table Tibble with cell_type/pos_markers/neg_markers.
#' @param marker_stats Output of fit_marker_stats().
#' @param cell_type_name Which cell type to plot.
#'
#' @return A ggraph/ggplot object (requires igraph + ggraph).
#' @export
plot_marker_priority_tree <- function(lineage_table, marker_stats, cell_type_name) {
  .assert_lineage_table(lineage_table)

  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("ggraph", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install igraph, ggraph, ggplot2 to use plot_marker_priority_tree().")
  }

  pos <- lineage_table$pos_markers[lineage_table$cell_type == cell_type_name][[1]]
  pos <- intersect(pos, names(marker_stats))
  if (length(pos) == 0) stop("No markers found for cell type: ", cell_type_name)

  weights <- vapply(pos, function(m) marker_stats[[m]]$weight %||% 0, numeric(1))
  ordered <- pos[order(weights, decreasing = TRUE)]

  df_edges <- data.frame(
    from = c("ROOT", ordered[-length(ordered)]),
    to   = ordered
  )
  df_nodes <- data.frame(name = c("ROOT", ordered))
  df_nodes$label <- df_nodes$name
  for (m in ordered) {
    w <- round(marker_stats[[m]]$weight %||% 0, 2)
    df_nodes$label[df_nodes$name == m] <- paste0(m, "\n(Sep: ", w, ")")
  }

  g <- igraph::graph_from_data_frame(df_edges, vertices = df_nodes)

  ggraph::ggraph(g, layout = "tree") +
    ggraph::geom_edge_link(arrow = grid::arrow(length = grid::unit(3, "mm")),
                           end_cap = ggraph::circle(3, "mm")) +
    ggraph::geom_node_point(size = 4) +
    ggraph::geom_node_label(ggplot2::aes(label = .data$label), repel = FALSE) +
    ggplot2::theme_void() +
    ggplot2::labs(title = paste("Soft Gating Priority Tree:", cell_type_name),
                  subtitle = "Markers ordered by separability weight")
}


#' Plot spatial probability map for a cell type
#'
#' Draws a spatial scatter plot colored by soft-gating probability
#' for a given cell type (e.g., "T cell", "Tumor").
#'
#' @param spe SpatialExperiment / SingleCellExperiment with P_* columns in colData.
#' @param cell_type Character. Must match the lineage_table$cell_type used in run_soft_gating().
#' @param image_col Column name for image ID (default "imageID").
#' @param x_col Coordinate column names (default auto-detect).
#' @param y_col Coordinate column names (default auto-detect).
#' @param point_size Size of points.
#'
#' @return A ggplot object.
#' @export
plot_probability_map <- function(spe,
                                 cell_type,
                                 image_col = "imageID",
                                 x_col = NULL,
                                 y_col = NULL,
                                 point_size = 0.8) {

  .assert_spe(spe)

  df <- as.data.frame(SummarizedExperiment::colData(spe))

  prob_col <- paste0("P_", gsub("\\s+", "_", cell_type))
  if (!prob_col %in% names(df)) {
    stop("Column ", prob_col, " not found in colData(spe). Did you run run_soft_gating()?")
  }

  coords <- tryCatch(SpatialExperiment::spatialCoords(spe), error = function(e) NULL)
  if (!is.null(coords)) {
    df$x <- coords[, 1]
    df$y <- coords[, 2]
    if (is.null(x_col)) x_col <- "x"
    if (is.null(y_col)) y_col <- "y"
  }

  if (is.null(x_col) || is.null(y_col) || !all(c(x_col, y_col) %in% names(df))) {
    stop("Could not determine x/y coordinates. Please specify x_col and y_col.")
  }

  if (!image_col %in% names(df)) {
    stop("image_col '", image_col, "' not found in colData(spe).")
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      color = .data[[prob_col]]
    )
  ) +
    ggplot2::geom_point(size = point_size, alpha = 0.7) +
    ggplot2::coord_fixed() +
    ggplot2::facet_wrap(stats::as.formula(paste("~", image_col))) +
    ggplot2::scale_color_viridis_c(option = "magma") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Probability map:", cell_type),
      color = "Probability",
      x = x_col,
      y = y_col
    )
}



#' Print a cell-type gating tree
#'
#' @param node A gating tree (node or leaf) produced by [build_fullcoverage_tree()].
#' @param indent String used internally for indentation during recursive printing.
#'
#' @return Invisibly returns `NULL`.
#' @export
print_celltype_tree <- function(node, indent = "") {
  if (node$type == "leaf") {
    cat(indent, "[leaf] depth =", node$depth,
        " n_cells =", length(node$cells), "\n", sep = "")
    return(invisible(NULL))
  }

  cat(indent,
      sprintf("[depth %d] %s > %.3f  (sep = %.2f, n = %d)\n",
              node$depth,
              node$marker,
              node$cutoff,
              node$sep_score %||% NA_real_,
              length(node$cells)),
      sep = "")

  # cat(indent, " ├─ low\n", sep = "")
  # print_celltype_tree(node$left, paste0(indent, " │  "))
  #
  # cat(indent, " └─ high\n", sep = "")
  # print_celltype_tree(node$right, paste0(indent, "    "))


  cat(indent, " |- low\n", sep = "")
  print_celltype_tree(node$left, paste0(indent, " |  "))

  cat(indent, " `- high\n", sep = "")
  print_celltype_tree(node$right, paste0(indent, "    "))
}

#' Convert a tree to a tibble representation
#'
#' @param node A gating tree (node or leaf).
#' @param parent Parent node id (used internally for recursion).
#' @param id Node id for the current subtree (used internally for recursion).
#'
#' @return A tibble with columns `id`, `parent`, `type`, `label`.
#' @export

tree_to_df <- function(node, parent = NA_character_, id = "root") {
  if (node$type == "leaf") {
    return(tibble::tibble(
      id = id,
      parent = parent,
      type = "leaf",
      label = paste0("leaf\nn=", length(node$cells))
    ))
  }

  this <- tibble::tibble(
    id = id,
    parent = parent,
    type = "node",
    label = sprintf("%s > %.2f\nsep=%.2f",
                    node$marker,
                    node$cutoff,
                    node$sep_score %||% NA_real_)
  )

  left <- tree_to_df(node$left, id, paste0(id, "_L"))
  right <- tree_to_df(node$right, id, paste0(id, "_R"))

  dplyr::bind_rows(this, left, right)
}

#' Plot a cell-type gating tree
#'
#' @param tree A gating tree produced by [build_fullcoverage_tree()].
#' @param title Plot title string.
#'
#' @return A `ggplot` object.
#' @export
plot_celltype_tree <- function(tree, title = "") {
  df <- tree_to_df(tree)

  g <- igraph::graph_from_data_frame(
    d = df |> dplyr::filter(!is.na(.data$parent)) |> dplyr::select(.data$parent, .data$id),
    vertices = df,
    directed = TRUE
  )

  ggraph::ggraph(g, layout = "dendrogram") +
    ggraph::geom_edge_elbow() +
    ggraph::geom_node_label(aes(label = .data$label),
                            size = 3,
                            label.size = 0.2,
                            fill = "white") +
    ggplot2::theme_void() +
    ggplot2::ggtitle(title)
}
