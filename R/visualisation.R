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

#' Plot marker intensity with optional GMM fit
#'
#' Visualizes a marker distribution, optionally with a two-component GMM fit
#' and cutoff overlays. Useful for QC and comparing marker behavior across
#' groups.
#'
#' @param spe A `SpatialExperiment`/`SingleCellExperiment`.
#' @param marker Marker name to plot (must match a row in the assay).
#' @param subtitle Subtitle for the plot.
#' @param binwidth Histogram bin width. If NULL, uses range/100.
#' @param label_col Optional column in `colData(spe)` used to subset to a
#'   specific cell type and to color the histogram.
#' @param cell_type Optional cell type value to subset when `label_col` is set.
#' @param do_fit Logical; fit a 2-component GMM and draw the components.
#' @param cutoff_method Passed to [fit_gmm_2()].
#' @param gmm_model_names Passed to [fit_gmm_2()].
#' @param label_levels Optional character vector of label levels to enforce for
#'   legend consistency.
#' @param label_colors Optional named or unnamed vector of colors to use for the
#'   label levels. If provided alongside `label_levels`, lengths must match.
#' @param drop_levels Logical; if FALSE, keep unused levels in the legend.
#'   Default FALSE.
#' @param assay_name Assay name to pull values from. Default "exprs".
#'
#' @return A `ggplot` object.
#' @export
plot_ct_marker_intensity <- function(spe,
                                     marker,
                                     subtitle,
                                     binwidth = NULL,
                                     label_col = NULL,
                                     cell_type = NULL,
                                     do_fit = TRUE,
                                     cutoff_method = "equal_posteriors",
                                     gmm_model_names = NULL,
                                     label_levels = NULL,
                                     label_colors = NULL,
                                     drop_levels = FALSE,
                                     assay_name = "exprs") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to use plot_ct_marker_intensity().")
  }

  .assert_spe(spe)

  if (!assay_name %in% SummarizedExperiment::assayNames(spe)) {
    stop("Assay '", assay_name, "' not found in spe.")
  }

  assay_mat <- SummarizedExperiment::assay(spe, assay_name)
  if (is.null(rownames(assay_mat)) || !marker %in% rownames(assay_mat)) {
    stop("Marker '", marker, "' not found in assay rownames.")
  }

  if (!is.null(label_col)) {
    if (!label_col %in% names(SummarizedExperiment::colData(spe))) {
      stop("label_col '", label_col, "' not found in colData(spe).")
    }
    if (!is.null(cell_type)) {
      spe <- spe[, SummarizedExperiment::colData(spe)[[label_col]] == cell_type]
      if (!ncol(spe)) stop("No cells found for the selected label subset.")
    }
  }

  marker_dist <- SummarizedExperiment::assay(spe, assay_name) |>
    t() |>
    as.data.frame() |>
    dplyr::pull(marker)

  fit <- NULL
  if (do_fit) {
    fit <- fit_gmm_2(
      marker_dist,
      cutoff_method = cutoff_method,
      gmm_model_names = gmm_model_names
    )
  }

  marker_df <- data.frame(marker = marker_dist, check.names = FALSE)

  resolve_binwidth <- function(x, binwidth) {
    if (!is.null(binwidth) && is.finite(binwidth) && binwidth > 0) {
      return(binwidth)
    }
    bw <- diff(range(x, na.rm = TRUE)) / 100
    if (!is.finite(bw) || bw <= 0) bw <- 0.01
    bw
  }
  if (!is.null(label_col)) {
    marker_df[[label_col]] <- as.character(SummarizedExperiment::colData(spe)[[label_col]])
    if (!is.null(label_levels)) {
      label_levels <- as.character(label_levels)
      marker_df[[label_col]] <- factor(marker_df[[label_col]], levels = label_levels)
    }

    label_levels_data <- marker_df[[label_col]]
    label_levels_data <- if (is.factor(label_levels_data)) levels(label_levels_data) else unique(as.character(label_levels_data))
    label_levels_use <- label_levels %||% label_levels_data

    if (is.null(label_colors)) {
      if (!requireNamespace("pals", quietly = TRUE)) {
        stop("Please install pals to use the default high-contrast palette in plot_ct_marker_intensity().")
      }
      label_colors <- as.vector(pals::polychrome(length(label_levels_use)))
      names(label_colors) <- label_levels_use
    } else {
      if (length(label_colors) != length(label_levels_use)) {
        stop("label_colors must be the same length as the label levels when provided.")
      }
      if (is.null(names(label_colors))) {
        names(label_colors) <- label_levels_use
      }
    }

    marker_df <- marker_df[is.finite(marker_df$marker) & !is.na(marker_df[[label_col]]), , drop = FALSE]
    if (!nrow(marker_df)) stop("No finite values found for marker.")

    bw <- resolve_binwidth(marker_df$marker, binwidth)

    p <- ggplot2::ggplot(marker_df) +
      ggplot2::aes(x = .data[["marker"]]) +
      ggplot2::geom_histogram(
        ggplot2::aes(fill = .data[[label_col]]),
        binwidth = bw,
        alpha = 0.6,
        position = "stack"
      ) +
      ggplot2::scale_fill_manual(values = label_colors, drop = drop_levels, limits = label_levels_use) +
      ggplot2::theme_classic()
  } else {
    marker_df <- marker_df[is.finite(marker_df$marker), , drop = FALSE]
    if (!nrow(marker_df)) stop("No finite values found for marker.")
    bw <- resolve_binwidth(marker_df$marker, binwidth)
    p <- ggplot2::ggplot(marker_df) +
      ggplot2::aes(x = .data[["marker"]]) +
      ggplot2::geom_histogram(
        binwidth = bw,
        fill = "grey70",
        color = "grey40",
        alpha = 1,
        show.legend = FALSE
      ) +
      ggplot2::theme_classic()
  }

  if (!is.null(fit)) {
    # cutoff <- fit$cutoff

    n_cells <- nrow(marker_df)
    bw <- resolve_binwidth(marker_df$marker, binwidth)

    x_seq <- seq(
      min(marker_df$marker, na.rm = TRUE),
      max(marker_df$marker, na.rm = TRUE),
      length.out = 200
    )

    gmm_df <- data.frame(
      x = x_seq,
      y1 = fit$p1 * dnorm(x_seq, mean = fit$mu1, sd = fit$s1),
      y2 = fit$p2 * dnorm(x_seq, mean = fit$mu2, sd = fit$s2)
    )
    gmm_df$y_total <- gmm_df$y1 + gmm_df$y2
    gmm_df$y1 <- gmm_df$y1 * n_cells * bw
    gmm_df$y2 <- gmm_df$y2 * n_cells * bw
    gmm_df$y_total <- gmm_df$y_total * n_cells * bw

    p <- p +
      ggplot2::geom_line(
        data = gmm_df,
        ggplot2::aes(x = .data[["x"]], y = .data[["y1"]]),
        linetype = "dashed",
        linewidth = 1,
        color = "red",
        show.legend = FALSE,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_line(
        data = gmm_df,
        ggplot2::aes(x = .data[["x"]], y = .data[["y2"]]),
        linetype = "dashed",
        linewidth = 1,
        color = "blue",
        show.legend = FALSE,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_line(
        data = gmm_df,
        ggplot2::aes(x = .data[["x"]], y = .data[["y_total"]]),
        linewidth = 0.5,
        alpha = 1,
        color = "black",
        show.legend = FALSE,
        inherit.aes = FALSE
      )
      # ggplot2::geom_vline(xintercept = cutoff, linetype = "solid", color = "orange", show.legend = FALSE) +
      # ggplot2::geom_vline(xintercept = fit$mu1, color = "red", show.legend = FALSE) +
      # ggplot2::geom_vline(xintercept = fit$mu2, color = "blue", show.legend = FALSE)
  }

  fill_label <- if (!is.null(label_col)) label_col else NULL

  p +
    ggplot2::labs(
      x = paste("Transformed and Normalised", marker, "expression"),
      y = "Count",
      title = paste(marker, "Histogram"),
      subtitle = paste(subtitle, "n =", ncol(spe)),
      fill = fill_label
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1))
}


#' Plot spatial probability map for a cell type
#'
#' Draws a spatial scatter plot colored by soft-gating probability
#' for a given cell type (e.g., "T cell", "Tumor").
#'
#' @param spe SpatialExperiment / SingleCellExperiment with P_* columns in colData.
#' @param cell_type Character. Must match the lineage_table$cell_type used in run_soft_gating().
#' @param image_col Column name for image ID (default "imageID").
#' @param image_index Numeric indices or character IDs of images to plot. Defaults
#'   to all images. If length > 1, results are faceted by image.
#' @param x_col Coordinate column names (default auto-detect).
#' @param y_col Coordinate column names (default auto-detect).
#' @param point_size Size of points.
#'
#' @return A ggplot object.
#' @export
plot_probability_map <- function(spe,
                                 cell_type,
                                 image_col = "imageID",
                                 image_index = NULL,
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

  image_values <- unique(stats::na.omit(df[[image_col]]))
  if (!length(image_values)) {
    stop("No image IDs found in column '", image_col, "'.")
  }

  if (is.null(image_index)) {
    selected_images <- image_values
  } else if (is.character(image_index)) {
    missing_imgs <- setdiff(image_index, image_values)
    if (length(missing_imgs)) {
      stop("image_index contains unknown image IDs: ", paste(missing_imgs, collapse = ", "))
    }
    selected_images <- image_index
  } else {
    if (!is.numeric(image_index) || any(is.na(image_index))) {
      stop("image_index must be numeric indices or character image IDs.")
    }
    image_index <- as.integer(image_index)
    if (any(image_index < 1L | image_index > length(image_values))) {
      stop("image_index values must be between 1 and ", length(image_values), ".")
    }
    selected_images <- image_values[image_index]
  }

  df <- df[df[[image_col]] %in% selected_images & !is.na(df[[image_col]]), , drop = FALSE]
  if (!nrow(df)) {
    stop("No cells found for the selected image(s).")
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      color = .data[[prob_col]]
    )
  ) +
    ggplot2::geom_point(size = point_size, alpha = 0.7) +
    ggplot2::coord_fixed() +
    ggplot2::scale_color_viridis_c(option = "magma") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Probability map:", cell_type),
      subtitle = if (length(unique(df[[image_col]])) > 1L) "Faceted by image" else paste0("Image: ", unique(df[[image_col]])),
      color = "Probability",
      x = x_col,
      y = y_col
    )

  if (length(unique(df[[image_col]])) > 1L) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", image_col)))
  }

  p
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


#' Plot a confusion matrix heatmap
#'
#' Creates a heatmap with counts overlaid for a confusion matrix (table or
#' matrix) of truth vs. predicted classes.
#'
#' @param conf_mat A confusion matrix as `table`, `matrix`, or data frame with
#'   columns `truth`, `pred`, and `Freq`.
#' @param title Plot title.
#' @param fill_low Low-end color for the fill gradient.
#' @param fill_high High-end color for the fill gradient.
#' @param subtitle Optional subtitle for the plot. If `NULL`, an automatic
#'   subtitle is generated (aggregation flag, CV folds if provided, and total `n`).
#' @param cv_folds Optional number of cross-validation folds to display in the subtitle.
#'   Use `NULL` to omit.
#' @param aggregated Logical. If `TRUE`, indicates the confusion matrix is summed/aggregated
#'   (used for subtitle text only).
#' @param plot_marginals Logical. If `TRUE`, add an extra "TOTAL" row/column showing marginal
#'   totals and the grand total.
#' @return A `ggplot` object.
#' @export
plot_confusion_matrix <- function(conf_mat,
                                  title = "Confusion matrix",
                                  subtitle = NULL,
                                  cv_folds = NULL,
                                  aggregated = TRUE,
                                  fill_low = "#f0f9e8",
                                  fill_high = "#08589e",
                                  plot_marginals = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to use plot_confusion_matrix().")
  }

  if (is.null(conf_mat)) {
    stop("conf_mat is NULL; no predictions available to plot.")
  }

  if (is.matrix(conf_mat) || is.table(conf_mat)) {
    df <- as.data.frame(as.table(conf_mat))
    names(df)[seq_len(3)] <- c("truth", "pred", "Freq")
  } else if (is.data.frame(conf_mat)) {
    needed <- c("truth", "pred", "Freq")
    missing_cols <- setdiff(needed, names(conf_mat))
    if (length(missing_cols)) {
      stop("conf_mat data.frame is missing columns: ", paste(missing_cols, collapse = ", "))
    }
    df <- conf_mat[, needed]
  } else {
    stop("conf_mat must be a table, matrix, or data.frame with truth/pred/Freq columns.")
  }

  if (!nrow(df)) stop("conf_mat is empty; nothing to plot.")
  if (all(df$Freq == 0, na.rm = TRUE)) stop("conf_mat has all zero counts; nothing to plot.")

  df$truth <- factor(df$truth)
  df$pred <- factor(df$pred, levels = levels(df$truth))

  # Prepare marginals data if requested
  marginals_df <- NULL
  if (plot_marginals) {
    # Calculate row totals (marginals for Truth)
    row_totals <- stats::aggregate(Freq ~ truth, data = df, FUN = sum)
    row_totals$pred <- "TOTAL"

    # Calculate column totals (marginals for Pred)
    col_totals <- stats::aggregate(Freq ~ pred, data = df, FUN = sum)
    col_totals$truth <- "TOTAL"

    # Grand total
    grand_total <- sum(df$Freq, na.rm = TRUE)
    grand_df <- data.frame(truth = "TOTAL", pred = "TOTAL", Freq = grand_total)

    # Combine marginals
    marginals_df <- rbind(row_totals, col_totals, grand_df)

    # Add factor levels for TOTAL
    df$truth <- factor(df$truth, levels = c(levels(df$truth), "TOTAL"))
    df$pred <- factor(df$pred, levels = c(levels(df$pred), "TOTAL"))
    marginals_df$truth <- factor(marginals_df$truth, levels = levels(df$truth))
    marginals_df$pred <- factor(marginals_df$pred, levels = levels(df$pred))
  }

  total_n <- sum(df$Freq, na.rm = TRUE)
  subtitle <- subtitle %||% {
    parts <- c(
      if (aggregated) "Summed confusion matrix" else NULL,
      if (!is.null(cv_folds)) paste0("CV folds = ", cv_folds) else NULL,
      paste0("n = ", total_n)
    )
    paste(parts, collapse = " | ")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$pred, y = .data$truth, fill = .data$Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$Freq), size = 3) +
    ggplot2::scale_fill_gradient(low = fill_low, high = fill_high) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Predicted", y = "Truth", fill = "Count") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

  # Add marginals as separate layer with fixed color
  if (plot_marginals && !is.null(marginals_df)) {
    p <- p +
      ggplot2::geom_tile(data = marginals_df, ggplot2::aes(x = .data$pred, y = .data$truth),
                         fill = "white", color = "gray50", linewidth = 1.2, show.legend = FALSE) +
      ggplot2::geom_text(data = marginals_df, ggplot2::aes(x = .data$pred, y = .data$truth, label = .data$Freq),
                         size = 3, color = "black", inherit.aes = FALSE)
  }

  p
}


#' Plot counts of a label column
#'
#' Horizontal bar plot of label counts with count annotations. Useful for
#' quickly checking class balance in `colData()` label columns.
#'
#' @param spe A `SpatialExperiment`/`SingleCellExperiment`.
#' @param label_col Column in `colData(spe)` containing labels to count.
#' @param title Plot title.
#' @param subtitle Optional subtitle. If no value provided, `label_col` is used.
#' @param include_na Logical; if TRUE, include missing values as "NA" level. Default FALSE.
#' @param show_legend Logical; show legend for bars. Default FALSE.
#'
#' @return A `ggplot` object.
#' @export
plot_label_counts <- function(spe,
                              label_col,
                              title = "Label counts",
                              subtitle = NULL,
                              include_na = FALSE,
                              show_legend = FALSE) {

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install ggplot2 and dplyr to use plot_label_counts().")
  }

  .assert_spe(spe)

  df <- SummarizedExperiment::colData(spe) |> as.data.frame()

  if (is.null(label_col) || !label_col %in% names(df)) {
    stop("label_col '", label_col, "' not found in colData(spe).")
  }

  labels <- df[[label_col]]
  if (!include_na) {
    keep <- !is.na(labels)
    labels <- labels[keep]
  }

  labels <- as.character(labels)
  labels[is.na(labels)] <- "NA"

  counts <- dplyr::tibble(label = labels) |>
    dplyr::count(.data$label, name = "n") |>
    dplyr::arrange(dplyr::desc(.data$n)) |>
    dplyr::mutate(label = factor(.data$label, levels = rev(.data$label)))

  subtitle <- subtitle %||% label_col

  ggplot2::ggplot(counts, ggplot2::aes(x = .data$label, y = .data$n, fill = .data$label)) +
    ggplot2::geom_col() +
    ggplot2::geom_label(ggplot2::aes(label = paste0("n = ", .data$n), y = .data$n),
                        hjust = -0.1, size = 3) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = "Count") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = if (show_legend) "right" else "none")
}


#' Plot agreement between two label columns
#'
#' Computes per-label agreement (Jaccard overlap) between two label columns in
#' `colData(spe)` and visualises it as a horizontal bar plot with class sizes.
#'
#' @param spe A `SpatialExperiment`/`SingleCellExperiment`.
#' @param label_col1 First label column name.
#' @param label_col2 Second label column name.
#' @param drop_na Logical; drop rows with NA in either label. Default TRUE.
#' @param title Plot title.
#' @param subtitle Optional subtitle.
#'
#' @return A `ggplot` object.
#' @export
plot_label_agreement <- function(spe,
                                 label_col1,
                                 label_col2,
                                 drop_na = TRUE,
                                 title = "Label agreement",
                                 subtitle = "total matches / union") {

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install ggplot2 and dplyr to use plot_label_agreement().")
  }

  agree_df <- label_agreement_rates(
    spe,
    label_col1 = label_col1,
    label_col2 = label_col2,
    drop_na = drop_na
  )

  if (!nrow(agree_df)) {
    stop("No labels available to compute agreement.")
  }

  agree_df <- agree_df |>
    dplyr::arrange(.data$agreement) |>
    dplyr::mutate(label = factor(.data$label, levels = .data$label))

  ggplot2::ggplot(agree_df, ggplot2::aes(x = .data$label, y = .data$agreement)) +
    ggplot2::geom_col(fill = "#3182bd") +
    # ggplot2::geom_text(ggplot2::aes(label = sprintf("n = %d", .data$union_n)),
    #                    hjust = -0.1, size = 3) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dotted", color = "red") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.25),
                                limits = c(0, 1.05),
                                expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(title = title,
                  subtitle = subtitle %||% paste0(label_col1, " vs ", label_col2),
                  x = NULL,
                  y = "Agreement (match / union)") +
    ggplot2::theme_classic()
}


#' Plot confusion matrix between two label columns
#'
#' Computes a confusion matrix for two label columns in `colData(spe)` and
#' visualises it with counts overlaid.
#'
#' @param spe A `SpatialExperiment`/`SingleCellExperiment`.
#' @param label_col1 First label column name.
#' @param label_col2 Second label column name.
#' @param drop_na Logical; drop rows with NA in either label. Default TRUE.
#' @param title Plot title.
#' @param subtitle Optional subtitle; defaults to "label_col1 vs label_col2".
#' @param fill_low Low-end color for the fill gradient. Default "#f0f9e8".
#' @param fill_high High-end color for the fill gradient. Default "#08589e".
#' @param plot_marginals Logical; if TRUE, plot marginal totals on the edges
#'   (right side and top) of the confusion matrix. Default FALSE.
#'
#' @return A `ggplot` object.
#' @export
plot_label_confusion_matrix <- function(spe,
                                        label_col1,
                                        label_col2,
                                        drop_na = TRUE,
                                        title = "Label confusion matrix",
                                        subtitle = NULL,
                                        fill_low = "#f0f9e8",
                                        fill_high = "#08589e",
                                        plot_marginals = FALSE) {

  .assert_spe(spe)

  conf <- label_confusion_matrix(
    spe,
    label_col1 = label_col1,
    label_col2 = label_col2,
    drop_na = drop_na
  )

  subtitle <- subtitle %||% paste0(label_col1, " vs ", label_col2)

  p <- plot_confusion_matrix(
    conf_mat = conf,
    title = title,
    subtitle = subtitle,
    aggregated = TRUE,
    fill_low = fill_low,
    fill_high = fill_high,
    plot_marginals = plot_marginals
  )

  p + ggplot2::labs(x = label_col2, y = label_col1)
}


#' Plot marker density across images
#'
#' Generates a density plot for a single marker using the specified assay, with
#' curves colored by the image/sample column. Useful for quick QC of
#' normalisation steps.
#'
#' @param spe A `SpatialExperiment`/`SingleCellExperiment` with the target assay
#'   and metadata in `colData()`.
#' @param marker Marker name to plot (must match a row in the assay).
#' @param assay_name Assay name to pull values from. Default "exprs".
#' @param image_col Column in `colData()` used to color densities. Default
#'   "image_name".
#' @param title Optional plot title. Defaults to "Density of \{marker\} (assay:
#'   \{assay_name\})".
#' @param show_legend Logical; whether to show the legend. Default `FALSE`.
#'
#' @return A `ggplot` object.
#' @export
plot_marker_density <- function(spe,
                                marker,
                                assay_name = "exprs",
                                image_col = "image_name",
                                title = NULL,
                                show_legend = FALSE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to use plot_marker_density().")
  }

  .assert_spe(spe)

  if (!assay_name %in% SummarizedExperiment::assayNames(spe)) {
    stop("Assay '", assay_name, "' not found in spe.")
  }

  assay_mat <- SummarizedExperiment::assay(spe, assay_name)

  if (is.null(rownames(assay_mat)) || !marker %in% rownames(assay_mat)) {
    stop("Marker '", marker, "' not found in assay '", assay_name, "'.")
  }

  meta <- SummarizedExperiment::colData(spe) |> as.data.frame()

  if (!image_col %in% names(meta)) {
    stop("image_col '", image_col, "' not found in colData(spe).")
  }

  meta[[marker]] <- as.numeric(assay_mat[marker, , drop = TRUE])

  title <- title %||% paste0("Density of ", marker, " (assay: ", assay_name, ")")

  ggplot2::ggplot(meta, ggplot2::aes(x = .data[[marker]], colour = .data[[image_col]])) +
    ggplot2::geom_density() +
    ggplot2::labs(
      title = title,
      x = marker,
      y = "Density",
      colour = image_col
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = if (show_legend) "right" else "none")
}


#' Plot per-class precision, recall, and F1
#'
#' Displays grouped bar plots for per-class precision, recall, and F1 with
#' support annotations. If a `fold` column is present (e.g., cross-validation
#' results), draws boxplots per class and facets by metric.
#'
#' @param class_metrics Data frame with columns `class`, `precision`, `recall`,
#'   `f1`, optionally `support`, and optionally `fold` for cross-validation.
#' @param title Plot title.
#' @param subtitle Optional subtitle for the plot. If `NULL`, an automatic
#' @return A `ggplot` object.
#' @export
plot_class_metrics <- function(class_metrics, title = "Per-class precision/recall/F1", subtitle = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE)) {
    stop("Please install ggplot2, dplyr, and tidyr to use plot_class_metrics().")
  }

  needed <- c("class", "precision", "recall", "f1")
  missing_cols <- setdiff(needed, names(class_metrics))
  if (length(missing_cols)) {
    stop("class_metrics is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  df <- class_metrics |>
    dplyr::mutate(class = factor(.data$class)) |>
    tidyr::pivot_longer(cols = c("precision", "recall", "f1"), names_to = "metric", values_to = "value")

  if (!nrow(df)) stop("class_metrics has no rows to plot.")

  has_folds <- "fold" %in% names(df)

  if (has_folds) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$class, y = .data$value, fill = .data$metric)) +
      ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.8), outlier.shape = 21, alpha = 0.6) +
      ggplot2::geom_jitter(ggplot2::aes(color = .data$metric),
                           position = ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                           alpha = 0.4, size = 1.5, show.legend = FALSE)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$class, y = .data$value, fill = .data$metric)) +
      ggplot2::geom_col(position = ggplot2::position_dodge())
  }

  if ("support" %in% names(class_metrics)) {
    supports <- class_metrics |>
      dplyr::select(.data$class, .data$support) |>
      dplyr::group_by(.data$class) |>
      dplyr::summarise(support = dplyr::first(.data$support), .groups = "drop") |>
      dplyr::mutate(class = factor(.data$class, levels = levels(df$class)))
    p <- p + ggplot2::geom_text(
      data = supports,
      ggplot2::aes(x = .data$class, y = 1.05, label = paste0("n = ", .data$support)),
      inherit.aes = FALSE,
      vjust = 0,
      size = 3
    )
  }

  total_n <- if ("support" %in% names(class_metrics)) sum(class_metrics$support, na.rm = TRUE) else NA_real_
  cv_txt <- if (has_folds) paste0("CV folds = ", length(unique(df$fold))) else NULL
  n_txt <- if (is.finite(total_n)) paste0("n = ", total_n) else NULL
  subtitle_auto <- paste(stats::na.omit(c(n_txt, cv_txt)), collapse = " | ")
  subtitle <- subtitle %||% subtitle_auto

  p <- p +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Class", y = "Score", fill = "Metric") +
    ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1.1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dotted", color = "red") +
    ggplot2::theme_classic() +
    ggplot2::coord_flip()

  p
}

#' Plot the cell type probabilities of a random sample of cells
#'
#' @param spe A spatial experiment object returned by `run_soft_gating()` or `run_tree_gating()` containing cell type probabilities in `colData()`. If omitted, tries `res$spe` when available.
#' @param image_index Index of the image from which to sample cells. Default 1.
#' @param sample_size Number of cells to sample (capped at available cells). Default 20.
#' @param image_col Column in `colData` holding the image IDs. Default "image_name".
#' @param flag_fn Optional function taking a numeric vector of probabilities for
#'   a cell and returning a single logical indicating confidence. Defaults to a
#'   Tukey-style rule matching [custom_labels()].
#' @param base_thresh Minimum threshold used by the default rule. Default 0.4.
#' @param quantile_cut Quantile used by the default rule. Default 0.75.
#' @param iqr_mult IQR multiplier for the default rule. Default 1.5.
#' @param summary_fun Function that returns a single numeric summary value per
#'   group (e.g., median or upper whisker). The value is plotted as y/ymin/ymax
#'   in `stat_summary`. Defaults to the same rule as the confidence flag
#'   (max(base_thresh, Q3 + iqr_mult * IQR)).
#'
#' @return a `ggplot` object
#' @export
plot_rand_cell_probs <- function(spe = NULL,
                                 image_index = 1,
                                 sample_size = 20,
                                 image_col = "image_name",
                                 flag_fn = NULL,
                                 base_thresh = 0.4,
                                 quantile_cut = 0.75,
                                 iqr_mult = 1.5,
                                 summary_fun = NULL) {

  # Fail fast if required packages are missing
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE)) {
    stop("Please install ggplot2, dplyr, and tidyr to use plot_rand_cell_probs().")
  }

  # Allow legacy default via res$spe while keeping an explicit argument path
  if (is.null(spe)) {
    res_obj <- get0("res", inherits = TRUE, ifnotfound = NULL)
    if (!is.null(res_obj) && !is.null(res_obj$spe)) {
      spe <- res_obj$spe
    } else {
      stop("Argument 'spe' is required when res$spe is not available.")
    }
  }

  .assert_spe(spe)

  df <- SummarizedExperiment::colData(spe) |> as.data.frame()

  if (!image_col %in% names(df)) {
    stop("image_col '", image_col, "' not found in colData(spe).")
  }

  prob_cols <- grep("^P_", names(df), value = TRUE)
  if (!length(prob_cols)) {
    stop("No probability columns starting with 'P_' found in colData(spe).")
  }

  image_values <- unique(stats::na.omit(df[[image_col]]))
  if (!length(image_values)) {
    stop("No image IDs found in column '", image_col, "'.")
  }
  if (image_index < 1 || image_index > length(image_values)) {
    stop("image_index must be between 1 and ", length(image_values), ".")
  }

  selected_image <- image_values[image_index]
  df <- df[df[[image_col]] == selected_image, , drop = FALSE]
  if (!nrow(df)) {
    stop("No cells found for image '", selected_image, "'.")
  }

  if (!"cell_number" %in% names(df)) {
    df$cell_number <- seq_len(nrow(df))
  }

  prob_mat <- as.matrix(df[, prob_cols, drop = FALSE])
  default_flag <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) return(FALSE)
    q <- stats::quantile(x, quantile_cut, na.rm = TRUE, names = FALSE)
    iqr <- stats::IQR(x, na.rm = TRUE)
    thr <- max(base_thresh, q + iqr_mult * iqr)
    max(x, na.rm = TRUE) > thr
  }
  default_summary_fun <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    q <- stats::quantile(x, quantile_cut, na.rm = TRUE, names = FALSE)
    iqr <- stats::IQR(x, na.rm = TRUE)
    max(base_thresh, q + iqr_mult * iqr)
  }
  flag_fn <- flag_fn %||% default_flag
  confident <- apply(prob_mat, 1, function(x) {
    out <- flag_fn(x)
    if (length(out) != 1 || !is.logical(out)) {
      stop("flag_fn must return a single logical value per cell.")
    }
    isTRUE(out)
  })
  df$confident <- confident

  if (sample_size < 1) {
    stop("sample_size must be >= 1.")
  }
  sample_size <- min(sample_size, nrow(df))

  sampled <- df[sample(seq_len(nrow(df)), sample_size), , drop = FALSE]

  summary_fun <- summary_fun %||% default_summary_fun
  if (!is.function(summary_fun)) {
    stop("summary_fun must be a function returning a single numeric value.")
  }

  sampled |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(prob_cols),
      names_to = "Cell_probs",
      names_prefix = "P_",
      values_to = "probability"
    ) |>
    dplyr::mutate(cell_number = factor(.data$cell_number)) |>
    ggplot2::ggplot(ggplot2::aes(x = .data$cell_number, y = .data$probability)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = NULL)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$Cell_probs, alpha = .data$confident), size = 1) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), guide = "none") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
    ggplot2::labs(title = "Cell type probabilities",
                  subtitle = paste0("Image: ", selected_image,
                                   " | Sampled cells: ", sample_size,
                                   " | Confident: ", sum(sampled$confident)),
                  x = "Cell",
                  y = "Probability") +
    ggplot2::stat_summary(
      fun.data = function(x) {
        val <- summary_fun(x)
        if (length(val) != 1 || !is.numeric(val) || is.na(val)) {
          stop("summary_fun must return a single non-NA numeric value.")
        }
        data.frame(y = val, ymin = val, ymax = val)
      },
      geom = "crossbar",
      color = "red",
      width = 0.5
    )
}

#' Plot cell-type probability histograms
#'
#' Faceted histograms of per-cell probabilities for each cell type, using the
#' probability columns already stored in `colData(spe)`
#' (e.g., `P_T_cell`, `P_Tumor`).
#'
#' @param res Result list containing an `spe` element with probability columns.
#' @param prob_prefix Prefix used for probability columns in `colData(spe)`. Default "P_".
#' @param binwidth Histogram bin width. Default 0.01.
#' @param drop_na Logical; drop NA probabilities before plotting. Default TRUE.
#' @param cutoff_fn Function taking a numeric vector of probabilities for a cell type and
#'   returning a single numeric cutoff to draw. Default [prob_quantile_cutoff()], which uses
#'   the 0.98 quantile.
#' @param label_col Optional column in `colData(spe)` used to color histograms
#'   by ground truth (or any label column).
#' @param label_levels Optional character vector of label levels to enforce for
#'   legend consistency.
#' @param label_colors Optional named or unnamed vector of colors to use for the
#'   label levels. If provided alongside `label_levels`, lengths must match.
#' @param drop_levels Logical; if FALSE, keep unused levels in the legend.
#'   Default FALSE.
#'
#' @return A `ggplot` object with faceted histograms.
#' @export
plot_score_hist <- function(res,
                           prob_prefix = "P_",
                           binwidth = 0.01,
                           drop_na = TRUE,
                           cutoff_fn = prob_quantile_cutoff,
                           label_col = NULL,
                           label_levels = NULL,
                           label_colors = NULL,
                           drop_levels = FALSE) {

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE)) {
    stop("Please install ggplot2, dplyr, and tidyr to use plot_score_hist().")
  }

  if (is.null(res) || is.null(res$spe)) {
    stop("Argument 'res' must contain an 'spe' element with probabilities in colData().")
  }

  spe <- res$spe

  .assert_spe(spe)

  meta <- SummarizedExperiment::colData(spe) |> as.data.frame()

  if (!"cell_id" %in% names(meta)) {
    meta$cell_id <- rownames(meta)
  }

  prob_cols <- names(meta)[startsWith(names(meta), prob_prefix)]
  if (!length(prob_cols)) {
    stop("No probability columns starting with '", prob_prefix, "' found in colData(spe).")
  }

  if (!is.null(label_col) && !label_col %in% names(meta)) {
    stop("label_col '", label_col, "' not found in colData(spe).")
  }

  plot_data <- meta[, unique(c("cell_id", prob_cols, label_col)), drop = FALSE]
  subtitle_txt <- paste0("All cells | n cells: ", length(unique(plot_data$cell_id)))

  prob_name_col <- if (!is.null(label_col) && label_col == "cell_type") {
    "prob_cell_type"
  } else {
    "cell_type"
  }

  plot_long <- plot_data |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(prob_cols),
      names_to = prob_name_col,
      names_prefix = prob_prefix,
      values_to = "scores"
    )

  if (drop_na) {
    plot_long <- dplyr::filter(plot_long, !is.na(.data$scores))
  }

  if (!is.null(label_col)) {
    plot_long <- plot_long[!is.na(plot_long[[label_col]]), , drop = FALSE]
  }

  if (!nrow(plot_long)) {
    stop("No probability values available to plot after filtering.")
  }

  if (!is.function(cutoff_fn)) {
    stop("cutoff_fn must be a function returning a single numeric cutoff per cell type.")
  }

  cutoff_df <- plot_long |>
    dplyr::group_by(.data[[prob_name_col]]) |>
    dplyr::summarise(
      cutoff = {
        val <- cutoff_fn(.data$scores)
        if (length(val) != 1 || !is.numeric(val)) {
          stop("cutoff_fn must return a single numeric value per cell type.")
        }
        as.numeric(val)
      },
      .groups = "drop"
    ) |>
    dplyr::filter(is.finite(.data$cutoff))

  use_label_col <- !is.null(label_col)

  if (use_label_col) {
    if (!is.null(label_levels)) {
      label_levels <- as.character(label_levels)
      plot_long[[label_col]] <- factor(plot_long[[label_col]], levels = label_levels)
    }

    label_levels_data <- plot_long[[label_col]]
    label_levels_data <- if (is.factor(label_levels_data)) levels(label_levels_data) else unique(as.character(label_levels_data))
    label_levels_use <- label_levels %||% label_levels_data

    if (is.null(label_colors)) {
      if (!requireNamespace("pals", quietly = TRUE)) {
        stop("Please install pals to use the default high-contrast palette in plot_score_hist().")
      }
      label_colors <- as.vector(pals::polychrome(length(label_levels_use)))
      names(label_colors) <- label_levels_use
    } else {
      if (length(label_colors) != length(label_levels_use)) {
        stop("label_colors must be the same length as the label levels when provided.")
      }
      if (is.null(names(label_colors))) {
        names(label_colors) <- label_levels_use
      }
    }
  }

  p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = .data$scores))

  if (use_label_col) {
    p <- p + ggplot2::geom_histogram(
      ggplot2::aes(fill = .data[[label_col]]),
      binwidth = binwidth,
      alpha = 0.6,
      boundary = 0,
      position = "stack"
    ) +
      ggplot2::scale_fill_manual(values = label_colors, drop = drop_levels, limits = label_levels_use)
  } else {
    p <- p + ggplot2::geom_histogram(binwidth = binwidth, alpha = 0.5, boundary = 0)
  }

  p +
    {if (nrow(cutoff_df)) ggplot2::geom_vline(data = cutoff_df, ggplot2::aes(xintercept = .data$cutoff),
                                              color = "red") else NULL} +
    ggplot2::facet_wrap(stats::as.formula(paste("~", prob_name_col)), scales = "free_y", axes = "all") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Histograms of Score for Each Cell Type",
                  subtitle = subtitle_txt,
            x = "Score",
                  y = "Count") +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

#' Barplot of positive label counts per cell
#'
#' Given a logical/boolean label matrix (cells x cell types), counts how many
#' labels are positive for each cell and plots the distribution from 0 up to
#' the maximum observed positives.
#'
#' @param label_mat Logical or numeric matrix with cells in rows and cell types
#'   in columns.
#' @param title Plot title. Default "Positive labels per cell".
#' @param xlab X-axis label. Default "Number of positive labels".
#' @param ylab Y-axis label. Default "Cell count".
#'
#' @return A `ggplot` object showing counts per cardinality.
#' @export
plot_label_cardinality <- function(label_mat,
                                   title = "Positive labels per cell",
                                   xlab = "Number of positive labels",
                                   ylab = "Cell count") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to use plot_label_cardinality().")
  }

  if (!is.matrix(label_mat)) stop("label_mat must be a matrix.")
  if (nrow(label_mat) == 0 || ncol(label_mat) == 0) {
    stop("label_mat must have at least one row and one column.")
  }

  # Convert to logical; treat non-finite as FALSE
  vals <- label_mat
  vals[!is.finite(vals)] <- 0
  vals <- vals != 0

  k <- rowSums(vals, na.rm = TRUE)
  if (!length(k)) stop("No rows available to summarize.")

  max_k <- max(k)
  k_vals <- 0:max_k
  counts <- tabulate(k + 1L, nbins = max_k + 1L)
  df <- data.frame(k = k_vals, n = counts)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$k, y = .data$n)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::scale_x_continuous(breaks = k_vals) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_classic()
}

#' Plots the spatial arrangement of classified cells
#'
#' Takes cell labels and coordinates from a SpatialExperiment; can plot a single
#' image or facet over multiple images with points colored by the provided
#' label column.
#'
#' @param spe SpatialExperiment object.
#' @param col_label Column in `colData(spe)` containing cell labels to color by.
#' @param image_index Index (or name) of the image(s) to view. Can be a single
#'   numeric, a numeric vector, or character vector of image IDs. Default 1.
#' @param image_col Column in `colData(spe)` containing image IDs. Default "imageID".
#' @param x_col Optional x-coordinate column in `colData(spe)`; auto-detected from
#'   `spatialCoords(spe)` when available.
#' @param y_col Optional y-coordinate column in `colData(spe)`; auto-detected from
#'   `spatialCoords(spe)` when available.
#' @param point_size Point size. Default 0.8.
#' @param label_levels Optional character vector of label levels to enforce for
#'   legend consistency. If provided, labels are coerced to a factor with these levels.
#' @param label_colors Optional named or unnamed vector of colors to use for the
#'   label levels. If provided alongside `label_levels`, lengths must match.
#' @param drop_levels Logical; if FALSE, keep unused levels in the legend.
#'   Default FALSE.
#' @param facet Logical; if TRUE, facet over selected images. If FALSE, plot a
#'   single image (first selected). Default TRUE.
#'
#' @return `ggplot` object
#' @export
plot_labelled_cells <- function(spe,
                                col_label,
                                image_index = 1,
                                image_col = "imageID",
                                x_col = NULL,
                                y_col = NULL,
                                point_size = 0.8,
                                label_levels = NULL,
                                label_colors = NULL,
                                drop_levels = FALSE,
                                facet = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to use plot_labelled_cells().")
  }

  .assert_spe(spe)

  df <- SummarizedExperiment::colData(spe) |> as.data.frame()

  if (is.null(col_label) || !col_label %in% names(df)) {
    stop("col_label '", col_label, "' not found in colData(spe).")
  }
  if (!image_col %in% names(df)) {
    stop("image_col '", image_col, "' not found in colData(spe).")
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

  image_values <- unique(stats::na.omit(df[[image_col]]))
  if (!length(image_values)) {
    stop("No image IDs found in column '", image_col, "'.")
  }

  if (is.character(image_index)) {
    missing_imgs <- setdiff(image_index, image_values)
    if (length(missing_imgs)) {
      stop("image_index contains unknown image IDs: ", paste(missing_imgs, collapse = ", "))
    }
    selected_images <- image_index
  } else {
    if (!is.numeric(image_index) || any(is.na(image_index))) {
      stop("image_index must be numeric indices or character image IDs.")
    }
    image_index <- as.integer(image_index)
    if (any(image_index < 1L | image_index > length(image_values))) {
      return(ggplot() + ggplot2::theme_void()) # return empty plot
      # stop("image_index values must be between 1 and ", length(image_values), ".")
    }
    selected_images <- image_values[image_index]
  }

  df_img <- df[df[[image_col]] %in% selected_images & !is.na(df[[image_col]]), , drop = FALSE]
  if (!nrow(df_img)) {
    stop("No cells found for the selected image(s).")
  }

  label_levels_data <- df_img[[col_label]]
  label_levels_data <- if (is.factor(label_levels_data)) levels(label_levels_data) else unique(as.character(label_levels_data))

  if (!is.null(label_levels)) {
    label_levels <- as.character(label_levels)
    df_img[[col_label]] <- factor(df_img[[col_label]], levels = label_levels)
  }

  label_levels_use <- label_levels %||% label_levels_data

  if (is.null(label_colors)) {
    if (!requireNamespace("pals", quietly = TRUE)) {
      stop("Please install pals to use the default high-contrast palette in plot_labelled_cells().")
    }
    label_colors <- as.vector(pals::polychrome(length(label_levels_use)))
    names(label_colors) <- label_levels_use
  } else {
    if (length(label_colors) != length(label_levels_use)) {
      stop("label_colors must be the same length as the label levels when provided.")
    }
    if (is.null(names(label_colors))) {
      names(label_colors) <- label_levels_use
    }
  }

  p <- ggplot2::ggplot(
    df_img,
    ggplot2::aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      color = .data[[col_label]]
    )
  ) +
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Spatial map of labelled cells",
      subtitle = if (facet) "Faceted by image" else paste0("Image: ", selected_images[1]),
      color = col_label,
      x = x_col,
      y = y_col
    )

  if (facet && length(selected_images) > 1L) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", image_col)))
  }

  p <- p + ggplot2::scale_color_manual(values = label_colors, drop = drop_levels, limits = label_levels_use) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))

  p
}


#' Plot label comparison dotplot
#'
#' Creates a dotplot comparing cell type distributions across multiple label
#' columns. The y-axis shows all unique cell types found across the specified
#' label columns (excluding Unknown/NA), x-axis shows counts, and points are
#' colored by which label column they came from.
#'
#' @param spe A `SpatialExperiment`/`SingleCellExperiment`.
#' @param label_cols Character vector of column names in `colData(spe)` to compare.
#'   Each column should contain cell type labels.
#' @param title Plot title.
#' @param subtitle Optional subtitle.
#' @param drop_unknown Logical; if TRUE, exclude "Unknown" from cell types. Default TRUE.
#' @param drop_na Logical; if TRUE, drop NA values. Default TRUE.
#' @param point_size Size of points. Default 4.
#' @param show_legend Logical; show legend. Default TRUE.
#'
#' @return A `ggplot` object.
#' @export
plot_label_dotplot <- function(spe,
                               label_cols,
                               title = "Label comparison",
                               subtitle = NULL,
                               drop_unknown = TRUE,
                               drop_na = TRUE,
                               point_size = 4,
                               show_legend = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE)) {
    stop("Please install ggplot2, dplyr, and tidyr to use plot_label_dotplot().")
  }

  .assert_spe(spe)

  if (length(label_cols) == 0) {
    stop("label_cols must contain at least one column name.")
  }

  df <- SummarizedExperiment::colData(spe) |> as.data.frame()

  missing_cols <- setdiff(label_cols, names(df))
  if (length(missing_cols)) {
    stop("Columns not found in colData(spe): ", paste(missing_cols, collapse = ", "))
  }

  # Gather counts from each label column
  count_list <- lapply(label_cols, function(col) {
    vals <- df[[col]]
    vals <- as.character(vals)

    if (drop_na) {
      vals <- vals[!is.na(vals)]
    }

    if (drop_unknown) {
      vals <- vals[vals != "Unknown"]
    }

    if (length(vals) == 0) return(NULL)

    counts <- table(vals)
    data.frame(
      cell_type = names(counts),
      total_count = as.numeric(counts),
      label_source = col,
      stringsAsFactors = FALSE
    )
  })

  count_df <- dplyr::bind_rows(count_list)

  if (!nrow(count_df)) {
    stop("No data available after filtering.")
  }

  # Order cell types by total counts across all sources
  cell_type_order <- count_df |>
    dplyr::group_by(.data$cell_type) |>
    dplyr::summarise(total = sum(.data$total_count), .groups = "drop") |>
    dplyr::arrange(.data$total) |>
    dplyr::pull(.data$cell_type)

  count_df$cell_type <- factor(count_df$cell_type, levels = cell_type_order)

  subtitle <- subtitle %||% paste0("Comparing ", length(label_cols), " label columns")

  ggplot2::ggplot(count_df, ggplot2::aes(x = .data$total_count,
                                        y = .data$cell_type,
                                        color = .data$label_source,
                                        shape = .data$label_source)) +
    ggplot2::geom_point(size = point_size, alpha = 0.5) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Total count (log10)",
      y = "Cell type",
      color = "Label source",
      shape = "Label source"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = if (show_legend) "right" else "none",
      axis.text.y = ggplot2::element_text(size = 10)
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_brewer(palette = "Set1", name = "Label source") +
    ggplot2::scale_shape_discrete(name = "Label source")
}

#' Plot pseudobulk marker heatmap
#'
#' Computes mean marker expression per label group and plots a heatmap
#' restricted to relevant markers defined in a lineage table.
#'
#' @param spe A `SpatialExperiment`/`SingleCellExperiment`.
#' @param label_col Column in `colData(spe)` containing labels.
#' @param lineage_table Tibble with cell_type/pos_markers/neg_markers.
#' @param assay_name Assay name to pull values from. Default "exprs".
#' @param cell_types Optional character vector to restrict label groups.
#' @param marker_groups Which marker sets to include from lineage_table.
#'   Default `c("pos_markers", "neg_markers")`.
#' @param include_all_markers Logical; if TRUE, use all assay markers rather than
#'   only lineage_table markers. Default TRUE.
#' @param drop_unknown Logical; drop "Unknown" labels. Default TRUE.
#' @param drop_na Logical; drop NA labels. Default TRUE.
#' @param cluster_rows Logical; if TRUE, cluster markers (rows). Default FALSE.
#' @param cluster_cols Logical; if TRUE, cluster labels (cols). Default FALSE.
#' @param scale Scaling for the heatmap. One of "none", "row", or "column".
#'   Default "none".
#' @param row_order Optional character vector of marker names to enforce when
#'   `cluster_rows = FALSE`. Any missing markers are ignored; any remaining markers
#'   not listed are appended in their existing order.
#' @param col_order Optional character vector of label names to enforce when
#'   `cluster_cols = FALSE`. Any missing labels are ignored; any remaining labels
#'   not listed are appended in their existing order.
#' @param highlight_lineage Logical; draw outlines for lineage markers per
#'   cell type. Positive markers are green, negative markers are red. Default FALSE.
#' @param heatmap_palette Palette name for the heatmap. Default "RdYlBu" for a
#'   blue-yellow-red diverging scale (pheatmap style). Set to a viridis option
#'   (e.g., "viridis", "magma") to use `ggplot2::scale_fill_viridis_c()`.
#' @param heatmap_direction Direction for viridis palettes. Default 1.
#' @param show_values Logical; draw mean values on tiles. Default TRUE.
#' @param value_digits Digits for tile labels. Default 2.
#' @param title Plot title. Default "Pseudobulk marker expression".
#' @param text_size Numeric. Text size for numeric values drawn on heatmap tiles
#'
#' @return A `ggplot` object.
#' @export
plot_pseudobulk_heatmap <- function(spe,
                                    label_col,
                                    lineage_table = NULL,
                                    assay_name = "exprs",
                                    cell_types = NULL,
                                    marker_groups = c("pos_markers", "neg_markers"),
                                    include_all_markers = TRUE,
                                    drop_unknown = TRUE,
                                    drop_na = TRUE,
                                    cluster_rows = FALSE,
                                    cluster_cols = FALSE,
                                    show_values = TRUE,
                                    value_digits = 2,
                                    title = "Pseudobulk marker expression",
                                    scale = "none",
                                    row_order = NULL,
                                    col_order = NULL,
                                    highlight_lineage = TRUE,
                                    heatmap_palette = "RdYlBu",
                                    heatmap_direction = 1,
                                    text_size = 2.5) {

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Please install ggplot2, dplyr, tidyr, and RColorBrewer to use plot_pseudobulk_heatmap().")
  }

  .assert_spe(spe)
  .assert_lineage_table(lineage_table)

  if (!assay_name %in% SummarizedExperiment::assayNames(spe)) {
    stop("Assay '", assay_name, "' not found in spe.")
  }

  meta <- SummarizedExperiment::colData(spe) |> as.data.frame()
  if (is.null(label_col) || !label_col %in% names(meta)) {
    stop("label_col '", label_col, "' not found in colData(spe).")
  }

  labels_raw <- meta[[label_col]]
  labels_raw <- as.character(labels_raw)

  keep <- rep(TRUE, length(labels_raw))
  if (drop_na) {
    keep <- keep & !is.na(labels_raw)
  }
  if (drop_unknown) {
    keep <- keep & labels_raw != "Unknown"
  }
  if (!is.null(cell_types)) {
    keep <- keep & labels_raw %in% cell_types
  }

  labels <- labels_raw[keep]

  if (!length(labels)) {
    stop("No labels available after filtering.")
  }

  missing_groups <- setdiff(marker_groups, c("pos_markers", "neg_markers"))
  if (length(missing_groups)) {
    stop("marker_groups must be any of: pos_markers, neg_markers.")
  }

  lineage_filtered <- lineage_table |>
    dplyr::filter(.data$cell_type %in% unique(labels))

  if (!nrow(lineage_filtered)) {
    stop("No matching cell types found in lineage_table for selected labels.")
  }

  assay_mat <- SummarizedExperiment::assay(spe, assay_name)
  if (is.null(rownames(assay_mat))) {
    stop("Assay has no rownames; cannot match markers.")
  }

  if (include_all_markers) {
    marker_list <- rownames(assay_mat)
  } else {
    marker_list <- lineage_filtered[, marker_groups, drop = FALSE] |>
      unlist(recursive = TRUE, use.names = FALSE)
    marker_list <- unique(marker_list)

    if (!length(marker_list)) {
      stop("No markers found in lineage_table for the selected marker_groups.")
    }

    marker_list <- intersect(marker_list, rownames(assay_mat))
    if (!length(marker_list)) {
      stop("No lineage markers found in assay rows for assay '", assay_name, "'.")
    }
  }

  assay_mat <- assay_mat[marker_list, keep, drop = FALSE]
  df <- t(assay_mat) |> as.data.frame()
  df[[label_col]] <- labels

  summary_df <- df |>
    dplyr::group_by(.data[[label_col]]) |>
    dplyr::summarise(
      dplyr::across(dplyr::where(is.numeric), mean, na.rm = TRUE),
      .groups = "drop"
    )

  # Determine ordering for rows/cols
  label_levels <- summary_df[[label_col]]
  marker_levels <- setdiff(names(summary_df), label_col)

  mat <- as.matrix(summary_df[, marker_levels, drop = FALSE])
  rownames(mat) <- summary_df[[label_col]]
  mat <- t(mat)
  rownames(mat) <- marker_levels
  colnames(mat) <- label_levels

  scale <- match.arg(scale, c("none", "row", "column"))

  apply_order <- function(order_vec, existing) {
    order_vec <- unique(order_vec)
    ordered <- intersect(order_vec, existing)
    remaining <- setdiff(existing, ordered)
    c(ordered, remaining)
  }

  if (cluster_rows && nrow(mat) > 1) {
    row_order <- rownames(mat)[stats::hclust(stats::dist(mat))$order]
  } else if (!is.null(row_order)) {
    if (!is.character(row_order)) {
      stop("row_order must be a character vector of marker names.")
    }
    row_order <- apply_order(row_order, rownames(mat))
  } else {
    row_order <- marker_levels
  }

  if (cluster_cols && ncol(mat) > 1) {
    col_order <- colnames(mat)[stats::hclust(stats::dist(t(mat)))$order]
  } else if (!is.null(col_order)) {
    if (!is.character(col_order)) {
      stop("col_order must be a character vector of label names.")
    }
    col_order <- apply_order(col_order, colnames(mat))
  } else {
    col_order <- label_levels
  }

  mat <- mat[row_order, col_order, drop = FALSE]

  zscore <- function(x) {
    mu <- mean(x, na.rm = TRUE)
    sigma <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sigma) || sigma == 0) {
      return(setNames(rep(0, length(x)), names(x)))
    }
    (x - mu) / sigma
  }

  if (scale != "none") {
    if (scale == "row") {
      mat <- t(apply(mat, 1, zscore))
    } else {
      mat <- apply(mat, 2, zscore)
    }
    mat <- as.matrix(mat)
    mat[is.na(mat)] <- 0
  }

  fill_label <- if (scale == "none") "Mean expr" else paste0("Scaled mean expr (", scale, ")")

  plot_df <- as.data.frame(mat) |>
    cbind(marker = rownames(mat)) |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of("marker"),
      names_to = "label",
      values_to = "expression"
    ) |>
    dplyr::mutate(
      marker = factor(.data$marker, levels = row_order),
      label = factor(.data$label, levels = col_order)
    )

  highlight_df <- NULL
  if (highlight_lineage) {
    extract_markers <- function(x) {
      unique(stats::na.omit(unlist(x, use.names = FALSE)))
    }

    pos_map <- lineage_filtered |>
      dplyr::select(.data$cell_type, .data$pos_markers)
    neg_map <- lineage_filtered |>
      dplyr::select(.data$cell_type, .data$neg_markers)

    pos_rows <- lapply(seq_len(nrow(pos_map)), function(i) {
      ct <- pos_map$cell_type[i]
      markers <- intersect(extract_markers(pos_map$pos_markers[i]), row_order)
      if (!length(markers)) return(NULL)
      data.frame(label = ct, marker = markers, outline_color = "#00ff00")
    })

    neg_rows <- lapply(seq_len(nrow(neg_map)), function(i) {
      ct <- neg_map$cell_type[i]
      markers <- intersect(extract_markers(neg_map$neg_markers[i]), row_order)
      if (!length(markers)) return(NULL)
      data.frame(label = ct, marker = markers, outline_color = "#ff00f7")
    })

    highlight_df <- dplyr::bind_rows(neg_rows, pos_rows)
    if (nrow(highlight_df)) {
      highlight_df$marker <- factor(highlight_df$marker, levels = row_order)
      highlight_df$label <- factor(highlight_df$label, levels = col_order)
    }
  }

  heatmap_colors <- grDevices::colorRampPalette(
    rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))
  )(100)

  value_range <- range(plot_df$expression, na.rm = TRUE, finite = TRUE)
  if (scale == "none") {
    fill_limits <- value_range
    midpoint <- stats::median(plot_df$expression, na.rm = TRUE)
  } else {
    max_abs <- max(abs(value_range))
    fill_limits <- c(-max_abs, max_abs)
    midpoint <- 0
  }

  use_viridis <- heatmap_palette %in% c(
    "viridis", "magma", "plasma", "inferno", "cividis", "rocket", "mako", "turbo"
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$label, y = .data$marker, fill = .data$expression)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    (if (use_viridis) {
      ggplot2::scale_fill_viridis_c(
        option = heatmap_palette,
        direction = heatmap_direction,
        na.value = "#f5f5f5"
      )
    } else {
      ggplot2::scale_fill_gradientn(
        colors = heatmap_colors,
        limits = fill_limits,
        na.value = "#f5f5f5"
      )
    }) +
    ggplot2::labs(title = title,
                  x = label_col,
                  y = "Marker",
                  fill = fill_label) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.height = ggplot2::unit(0.8, "cm")
    )

  if (!is.null(highlight_df) && nrow(highlight_df)) {
    p <- p + ggplot2::geom_tile(
      data = highlight_df,
      ggplot2::aes(x = .data$label, y = .data$marker),
      color = highlight_df$outline_color,
      fill = NA,
      linewidth = 0.8,
      inherit.aes = FALSE
    )
  }

  if (show_values) {
    span <- max(abs(value_range - midpoint))
    if (!is.finite(span) || span == 0) {
      plot_df$text_color <- "black"
    } else {
      # Use a mid-range band for black text; dark extremes use white for contrast.
      plot_df$text_color <- ifelse(abs(plot_df$expression - midpoint) <= 100 * span, "black", "white")
    }

    p <- p + ggplot2::geom_text(
      data = plot_df,
      ggplot2::aes(
        label = round(.data$expression, value_digits)
      ),
      color = plot_df$text_color,
      size = text_size,
      fontface = "bold",
      show.legend = FALSE
    )
  }
  print(p)

  return(
    list(
      plot = p, plot_df = plot_df
    )
  )
}
