## Create the cytoGateR_example package dataset
##
## Run this script from the cytoGateR package source directory after installing
## cytoGateR and its suggested data dependency, imcdatasets:
##
##   Rscript inst/scripts/create_cytoGateR_example.R
##
## An alternative output path can be supplied as the first command-line
## argument. The default is data/cytoGateR_example.rda.

required_packages <- c(
  "BiocParallel",
  "cytoGateR",
  "imcdatasets",
  "SingleCellExperiment",
  "SpatialExperiment",
  "SummarizedExperiment",
  "S4Vectors",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the packages required to create the example data: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

arguments <- commandArgs(trailingOnly = TRUE)
output_file <- if (length(arguments)) {
  arguments[[1L]]
} else {
  file.path("data", "cytoGateR_example.rda")
}

if (!length(arguments) &&
    (!file.exists("DESCRIPTION") || !dir.exists("data"))) {
  stop(
    "Run this script from the cytoGateR package source directory, or supply ",
    "an explicit output path as the first command-line argument.",
    call. = FALSE
  )
}

## Obtain the source SpatialExperiment from the imcdatasets ExperimentHub
## resource. The download is cached by the Bioconductor infrastructure.
spe <- imcdatasets::IMMUcan_2022_CancerExample("spe")

## Use two complete samples to retain realistic spatial and biological
## variation while keeping the installed dataset compact.
selected_samples <- c("Patient4_006", "Patient4_008")
subspe <- spe[, spe$sample_id %in% selected_samples]

## Retain only information used by package examples and vignettes.
SummarizedExperiment::assays(subspe) <-
  SummarizedExperiment::assays(subspe)["exprs"]
SingleCellExperiment::reducedDims(subspe) <- S4Vectors::SimpleList()
S4Vectors::metadata(subspe) <- list()
SummarizedExperiment::rowData(subspe) <-
  S4Vectors::DataFrame(row.names = rownames(subspe))

lineage_table <- tibble::tibble(
  cell_type = c(
    "Tumor", "Stroma", "Myeloid", "Neutrophil", "Bcell",
    "Plasma_cell", "BnTcell", "CD8", "CD4", "Treg"
  ),
  pos_markers = list(
    c("CDH1", "CA9", "PD_L1"),
    c("SMA", "PDGFRB"),
    c(
      "CD68", "CD163", "CD14", "CD11c", "CD33", "CD206", "CD303",
      "HLA_DR", "IDO1", "VISTA", "CD40"
    ),
    c("MPO", "CD15", "CD16"),
    c("CD20", "CD40", "HLA_DR", "B2M"),
    c("CD38", "CD27"),
    c("CD3e", "CD20", "CD7", "B2M"),
    c("CD3e", "CD8a", "GZMB", "TCF7"),
    c("CD3e", "CD4"),
    c("CD3e", "CD4", "FOXP3", "ICOS")
  ),
  neg_markers = list(
    c("CD3e", "CD20", "SMA", "CD68"),
    c("CDH1", "CD3e", "CD20", "MPO"),
    c("CD3e", "CD20", "CDH1"),
    c("CD3e", "CD20"),
    c("CD3e", "MPO", "CDH1"),
    c("CD20", "CD3e", "CD4"),
    c("CDH1", "SMA", "MPO"),
    c("CD4", "FOXP3", "CD20"),
    c("CD8a", "FOXP3", "CD20"),
    c("CD8a", "CD20")
  )
)

## Tree gating is performed once while preparing the package data rather than
## during examples or vignette construction.
set.seed(123)
tree_result <- cytoGateR::run_tree_gating(
  subspe,
  lineage_table,
  assay_name = "exprs",
  max_depth = 4,
  min_cells = 100,
  min_score = 0.05,
  uncert_thresh = 0.25,
  neg_strength = 0.5,
  cutoff_method = "equal_posteriors",
  workers = 1
)

tree_result <- cytoGateR::apply_cutoff_labels(
  tree_result,
  cutoff_fn = function(x) {
    cytoGateR::prob_quantile_cutoff(x, prob = 0.95)
  },
  label_col = "core_group",
  unknown_label = "Unassigned"
)

features <- unique(unlist(
  lineage_table$pos_markers,
  use.names = FALSE
))
features <- intersect(features, rownames(tree_result$spe))

## Train and apply a compact serial weighted-kNN model.
knn_reference <- cytoGateR::train_custom_knn(
  spe = tree_result$spe,
  assay_name = "exprs",
  label_col = "core_group",
  unknown_label = "Unassigned",
  features = features,
  cv_folds = 3,
  repeats = 2,
  agreement_thresh = 0.6,
  k = 3,
  method = "pearson",
  seed = 123,
  chunk_size = 500L,
  BPPARAM = BiocParallel::SerialParam()
)

knn_result <- cytoGateR::predict_unknown_with_knn(
  spe = knn_reference$spe,
  knn_ref = knn_reference,
  assay_name = "exprs",
  label_col = "cleaned_core_label",
  out_col = "knn_label",
  pred_col = "knn_pred",
  unknown_label = "Unassigned",
  unassigned_label = "Unassigned",
  threshold = 0.5,
  k = 3,
  dist_method = "pearson",
  chunk_size = 500L
)
subspe <- knn_result$spe

## Combine probabilities for retained reference cells and predicted cells.
full_knn_prob_mat <- rbind(
  knn_reference$core_prob_mat,
  knn_result$prob_mat
)

if (anyDuplicated(rownames(full_knn_prob_mat))) {
  stop("Duplicated cell names found while combining kNN probabilities.")
}

missing_cells <- setdiff(colnames(subspe), rownames(full_knn_prob_mat))
extra_cells <- setdiff(rownames(full_knn_prob_mat), colnames(subspe))
if (length(missing_cells) || length(extra_cells)) {
  stop(
    "Combined kNN probabilities do not match the cells in the example: ",
    length(missing_cells), " missing and ",
    length(extra_cells), " extra."
  )
}

full_knn_prob_mat <-
  full_knn_prob_mat[colnames(subspe), , drop = FALSE]

if (anyNA(full_knn_prob_mat) ||
    !isTRUE(all.equal(
      unname(rowSums(full_knn_prob_mat)),
      rep(1, nrow(full_knn_prob_mat)),
      tolerance = 1e-8
    ))) {
  stop("The complete kNN probability matrix is not normalized.")
}

for (cell_type_name in colnames(full_knn_prob_mat)) {
  output_name <- paste0(
    "KNN_P_",
    gsub("\\s+", "_", cell_type_name)
  )
  SummarizedExperiment::colData(subspe)[[output_name]] <-
    full_knn_prob_mat[, cell_type_name]
}

tree_probability_columns <- grep(
  "^P_",
  names(SummarizedExperiment::colData(subspe)),
  value = TRUE
)
knn_probability_columns <- grep(
  "^KNN_P_",
  names(SummarizedExperiment::colData(subspe)),
  value = TRUE
)

keep_coldata <- c(
  "sample_id",
  "image_name",
  "cell_type",
  "cell_type_hard",
  "prob_best",
  "core_group",
  "knn_label",
  "knn_confidence",
  tree_probability_columns,
  knn_probability_columns
)
SummarizedExperiment::colData(subspe) <-
  SummarizedExperiment::colData(subspe)[, keep_coldata, drop = FALSE]

## Remove unused factor levels inherited from samples not retained.
for (column_name in names(SummarizedExperiment::colData(subspe))) {
  column_value <- SummarizedExperiment::colData(subspe)[[column_name]]
  if (is.factor(column_value)) {
    SummarizedExperiment::colData(subspe)[[column_name]] <-
      droplevels(column_value)
  }
}

cytoGateR_example <- subspe
output_directory <- dirname(output_file)
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

save(
  cytoGateR_example,
  file = output_file,
  compress = "xz",
  compression_level = 9
)

message(
  "Saved ", ncol(cytoGateR_example), " cells and ",
  nrow(cytoGateR_example), " markers to ", normalizePath(output_file)
)
