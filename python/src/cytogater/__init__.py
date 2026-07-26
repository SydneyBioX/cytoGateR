"""Marker-aware cell annotation for Python."""

from .core import (
    assign_soft_labels,
    build_fullcoverage_tree,
    calculate_soft_scores,
    collect_path_scores,
    fit_gmm_2,
    fit_marker_stats,
    gmm_equal_posterior_cutoff,
    marker_separability,
    neg_penalty,
    run_soft_gating,
    run_tree_gating,
    score_marker_logistic,
    score_marker_rank,
    tree_prob,
)
from .gating.hierarchical import (
    build_hierarchical_reference,
    build_lineage_hierarchy,
)
from .metrics import (
    apply_cutoff_labels,
    assign_confident_labels,
    calculate_f1,
    class_metrics_from_fit,
    compute_custom_labels,
    custom_labels,
    label_agreement_rates,
    label_confusion_matrix,
    prob_mad_cutoff,
    prob_quantile_cutoff,
    probability_label_matrix,
)
from .models.neural_network import predict_unknown_with_dl, train_custom_dl
from .plotting import (
    plot_celltype_tree,
    plot_class_metrics,
    plot_confusion_matrix,
    plot_ct_marker_intensity,
    plot_label_agreement,
    plot_label_cardinality,
    plot_label_confusion_matrix,
    plot_label_counts,
    plot_label_dotplot,
    plot_labelled_cells,
    plot_marker_density,
    plot_marker_priority_tree,
    plot_probability_map,
    plot_pseudobulk_heatmap,
    plot_rand_cell_probs,
    plot_score_hist,
    print_celltype_tree,
    tree_to_df,
)

try:
    from .models.knn import (
        predict_hierarchical_knn_recursive,
        predict_unknown_with_knn,
        predict_wknn_multi,
        train_custom_knn,
    )
    from .models.random_forest import (
        predict_unknown_with_randomforest,
        rf_metric_table,
        rf_metric_text,
        train_custom_randomforest,
    )
    from .uncertainty import calculate_spatial_prior_labels, calculate_uncertainty
except ImportError:
    # Gating remains importable in lightweight environments. A normal package
    # installation installs scikit-learn and exposes these functions.
    pass

__all__ = [
    "assign_soft_labels",
    "build_fullcoverage_tree",
    "calculate_soft_scores",
    "collect_path_scores",
    "fit_gmm_2",
    "fit_marker_stats",
    "gmm_equal_posterior_cutoff",
    "marker_separability",
    "neg_penalty",
    "run_soft_gating",
    "run_tree_gating",
    "score_marker_logistic",
    "score_marker_rank",
    "tree_prob",
    "apply_cutoff_labels", "assign_confident_labels",
    "build_hierarchical_reference", "build_lineage_hierarchy",
    "calculate_f1", "calculate_spatial_prior_labels", "calculate_uncertainty",
    "class_metrics_from_fit", "compute_custom_labels", "custom_labels",
    "label_agreement_rates", "label_confusion_matrix",
    "predict_hierarchical_knn_recursive", "predict_unknown_with_dl",
    "predict_unknown_with_knn", "predict_unknown_with_randomforest",
    "predict_wknn_multi", "prob_mad_cutoff", "prob_quantile_cutoff",
    "probability_label_matrix", "rf_metric_table", "rf_metric_text",
    "train_custom_dl", "train_custom_knn", "train_custom_randomforest",
    "plot_celltype_tree", "plot_class_metrics", "plot_confusion_matrix",
    "plot_ct_marker_intensity", "plot_label_agreement",
    "plot_label_cardinality", "plot_label_confusion_matrix",
    "plot_label_counts", "plot_label_dotplot", "plot_labelled_cells",
    "plot_marker_density", "plot_marker_priority_tree",
    "plot_probability_map", "plot_pseudobulk_heatmap",
    "plot_rand_cell_probs", "plot_score_hist", "print_celltype_tree",
    "tree_to_df",
]
