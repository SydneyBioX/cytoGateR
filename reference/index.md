# Package index

## Marker-aware gating

Build interpretable cell-type scores, gates, and gating trees.

- [`fit_gmm_2()`](https://sydneybiox.github.io/cytoGateR/reference/fit_gmm_2.md)
  : Fit a 2-component Gaussian mixture model
- [`fit_marker_stats()`](https://sydneybiox.github.io/cytoGateR/reference/fit_marker_stats.md)
  : Fit per-marker 2-GMM statistics for soft gating
- [`gmm_equal_posterior_cutoff()`](https://sydneybiox.github.io/cytoGateR/reference/gmm_equal_posterior_cutoff.md)
  : Compute the cutoff where posteriors are equal
- [`marker_separability()`](https://sydneybiox.github.io/cytoGateR/reference/marker_separability.md)
  : Marker separability
- [`neg_penalty()`](https://sydneybiox.github.io/cytoGateR/reference/neg_penalty.md)
  : Negative marker penalty
- [`score_marker_logistic()`](https://sydneybiox.github.io/cytoGateR/reference/score_marker_logistic.md)
  : Logistic marker score
- [`score_marker_rank()`](https://sydneybiox.github.io/cytoGateR/reference/score_marker_rank.md)
  : Rank based marker score
- [`calculate_soft_scores()`](https://sydneybiox.github.io/cytoGateR/reference/calculate_soft_scores.md)
  : Calculate soft gating probabilities per cell type
- [`assign_soft_labels()`](https://sydneybiox.github.io/cytoGateR/reference/assign_soft_labels.md)
  : Assign hard labels from probability matrix
- [`run_soft_gating()`](https://sydneybiox.github.io/cytoGateR/reference/run_soft_gating.md)
  : Run soft-gating pipeline end-to-end
- [`run_tree_gating()`](https://sydneybiox.github.io/cytoGateR/reference/run_tree_gating.md)
  : Run tree-based soft gating on a SpatialExperiment
- [`build_fullcoverage_tree()`](https://sydneybiox.github.io/cytoGateR/reference/build_fullcoverage_tree.md)
  : Build full-coverage gating tree
- [`collect_path_scores()`](https://sydneybiox.github.io/cytoGateR/reference/collect_path_scores.md)
  : Collect scores along a tree path
- [`tree_prob()`](https://sydneybiox.github.io/cytoGateR/reference/tree_prob.md)
  : Tree probability for one cell

## Reference cells and prediction

Select confident references and classify cells that remain unresolved.

- [`apply_cutoff_labels()`](https://sydneybiox.github.io/cytoGateR/reference/apply_cutoff_labels.md)
  : Apply cutoff-based labels to a result list
- [`assign_confident_labels()`](https://sydneybiox.github.io/cytoGateR/reference/assign_confident_labels.md)
  : Assign labels with uncertainty handling
- [`compute_custom_labels()`](https://sydneybiox.github.io/cytoGateR/reference/compute_custom_labels.md)
  : Flag confident cells from a probability matrix
- [`custom_labels()`](https://sydneybiox.github.io/cytoGateR/reference/custom_labels.md)
  : Add custom labels to a SpatialExperiment
- [`prob_mad_cutoff()`](https://sydneybiox.github.io/cytoGateR/reference/prob_mad_cutoff.md)
  : Median + MAD cutoff helper
- [`prob_quantile_cutoff()`](https://sydneybiox.github.io/cytoGateR/reference/prob_quantile_cutoff.md)
  : Default quantile-based cutoff
- [`train_custom_knn()`](https://sydneybiox.github.io/cytoGateR/reference/train_custom_knn.md)
  : Train kNN Reference for Cell Type Prediction
- [`predict_unknown_with_knn()`](https://sydneybiox.github.io/cytoGateR/reference/predict_unknown_with_knn.md)
  : Predict and fill unknown labels using a weighted kNN model
- [`predict_wknn_multi()`](https://sydneybiox.github.io/cytoGateR/reference/predict_wknn_multi.md)
  : Multi-Metric Weighted kNN
- [`train_custom_randomforest()`](https://sydneybiox.github.io/cytoGateR/reference/train_custom_randomforest.md)
  : Train a custom random forest with consensus label cleaning
- [`predict_unknown_with_randomforest()`](https://sydneybiox.github.io/cytoGateR/reference/predict_unknown_with_randomforest.md)
  : Predict and fill unknown labels using a trained random forest model
- [`train_custom_dl()`](https://sydneybiox.github.io/cytoGateR/reference/train_custom_dl.md)
  : Train a custom deep-learning classifier with consensus label
  cleaning
- [`predict_unknown_with_dl()`](https://sydneybiox.github.io/cytoGateR/reference/predict_unknown_with_dl.md)
  : Predict unknown cells using a trained deep-learning model
- [`build_lineage_hierarchy()`](https://sydneybiox.github.io/cytoGateR/reference/build_lineage_hierarchy.md)
  : Build Lineage Hierarchy from Core Cells
- [`build_hierarchical_reference()`](https://sydneybiox.github.io/cytoGateR/reference/build_hierarchical_reference.md)
  : Build Hierarchical Reference with Flexible Marker Selection
- [`predict_hierarchical_knn_recursive()`](https://sydneybiox.github.io/cytoGateR/reference/predict_hierarchical_knn_recursive.md)
  : Recursive Hierarchical kNN with BiocParallel Ensemble

## Metrics and uncertainty

Quantify agreement, classification performance, and spatial uncertainty.

- [`calculate_f1()`](https://sydneybiox.github.io/cytoGateR/reference/calculate_f1.md)
  : Calculate F1 Score per Cell Type
- [`calculate_uncertainty()`](https://sydneybiox.github.io/cytoGateR/reference/calculate_uncertainty.md)
  : Calculate Full Uncertainty Suite with Sample Grouping
- [`calculate_spatial_prior_labels()`](https://sydneybiox.github.io/cytoGateR/reference/calculate_spatial_prior_labels.md)
  : Apply Protected Spatial Priors to kNN Probabilities
- [`class_metrics_from_fit()`](https://sydneybiox.github.io/cytoGateR/reference/class_metrics_from_fit.md)
  : Derive per-class metrics from a fitted model
- [`label_agreement_rates()`](https://sydneybiox.github.io/cytoGateR/reference/label_agreement_rates.md)
  : Agreement rates between two label columns
- [`label_confusion_matrix()`](https://sydneybiox.github.io/cytoGateR/reference/label_confusion_matrix.md)
  : Confusion matrix between two label columns
- [`probability_label_matrix()`](https://sydneybiox.github.io/cytoGateR/reference/probability_label_matrix.md)
  : Build a binary label matrix from probabilities
- [`rf_metric_table()`](https://sydneybiox.github.io/cytoGateR/reference/rf_metric_table.md)
  : Tabular and text summaries of RF performance
- [`rf_metric_text()`](https://sydneybiox.github.io/cytoGateR/reference/rf_metric_text.md)
  : Text summary of random forest performance

## Visualisation

Inspect marker distributions, spatial labels, probabilities, and trees.

- [`plot_celltype_tree()`](https://sydneybiox.github.io/cytoGateR/reference/plot_celltype_tree.md)
  : Plot a cell-type gating tree
- [`plot_class_metrics()`](https://sydneybiox.github.io/cytoGateR/reference/plot_class_metrics.md)
  : Plot per-class precision, recall, and F1
- [`plot_confusion_matrix()`](https://sydneybiox.github.io/cytoGateR/reference/plot_confusion_matrix.md)
  : Plot a confusion matrix heatmap
- [`plot_ct_marker_intensity()`](https://sydneybiox.github.io/cytoGateR/reference/plot_ct_marker_intensity.md)
  : Plot marker intensity with optional GMM fit
- [`plot_label_agreement()`](https://sydneybiox.github.io/cytoGateR/reference/plot_label_agreement.md)
  : Plot agreement between two label columns
- [`plot_label_cardinality()`](https://sydneybiox.github.io/cytoGateR/reference/plot_label_cardinality.md)
  : Barplot of positive label counts per cell
- [`plot_label_confusion_matrix()`](https://sydneybiox.github.io/cytoGateR/reference/plot_label_confusion_matrix.md)
  : Plot confusion matrix between two label columns
- [`plot_label_counts()`](https://sydneybiox.github.io/cytoGateR/reference/plot_label_counts.md)
  : Plot counts of a label column
- [`plot_label_dotplot()`](https://sydneybiox.github.io/cytoGateR/reference/plot_label_dotplot.md)
  : Plot label comparison dotplot
- [`plot_labelled_cells()`](https://sydneybiox.github.io/cytoGateR/reference/plot_labelled_cells.md)
  : Plots the spatial arrangement of classified cells
- [`plot_marker_density()`](https://sydneybiox.github.io/cytoGateR/reference/plot_marker_density.md)
  : Plot marker density across images
- [`plot_marker_priority_tree()`](https://sydneybiox.github.io/cytoGateR/reference/plot_marker_priority_tree.md)
  : Plot a marker priority "tree" for a given cell type (weights order)
- [`plot_probability_map()`](https://sydneybiox.github.io/cytoGateR/reference/plot_probability_map.md)
  : Plot spatial probability map for a cell type
- [`plot_pseudobulk_heatmap()`](https://sydneybiox.github.io/cytoGateR/reference/plot_pseudobulk_heatmap.md)
  : Plot pseudobulk marker heatmap
- [`plot_rand_cell_probs()`](https://sydneybiox.github.io/cytoGateR/reference/plot_rand_cell_probs.md)
  : Plot the cell type probabilities of a random sample of cells
- [`plot_score_hist()`](https://sydneybiox.github.io/cytoGateR/reference/plot_score_hist.md)
  : Plot cell-type probability histograms
- [`print_celltype_tree()`](https://sydneybiox.github.io/cytoGateR/reference/print_celltype_tree.md)
  : Print a cell-type gating tree
- [`tree_to_df()`](https://sydneybiox.github.io/cytoGateR/reference/tree_to_df.md)
  : Convert a tree to a tibble representation

## Data

- [`cytoGateR_example`](https://sydneybiox.github.io/cytoGateR/reference/cytoGateR_example.md)
  : Minimal spatial proteomics example
