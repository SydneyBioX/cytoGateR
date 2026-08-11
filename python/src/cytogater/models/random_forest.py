"""Random-forest consensus cleaning and unknown-cell prediction."""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold

from ..preprocessing import expression_frame


def train_custom_randomforest(
    spe,
    label_col="cutoff_label",
    assay_name="exprs",
    unknown_label="Unknown",
    num_trees=200,
    mtry=None,
    seed=None,
    features="all",
    cv_folds=5,
    repeats=10,
    agreement_thresh=0.8,
    num_threads=2,
    **kwargs,
):
    """Train the R-equivalent probability random forest.

    ``num.trees`` may be supplied through ``kwargs`` for direct R-call migration.
    """
    num_trees = kwargs.pop("num.trees", num_trees)
    frame = expression_frame(spe, assay_name)
    features_used = (
        list(frame.columns)
        if isinstance(features, str) and features == "all"
        else [x for x in features if x in frame]
    )
    labels = spe.obs[label_col]
    core = labels.notna() & (labels != unknown_label)
    x, y = frame.loc[core, features_used], labels[core].astype(str)
    class_counts = y.value_counts()
    eligible = y.isin(class_counts[class_counts >= 2].index).to_numpy()
    if not eligible.any():
        raise ValueError(
            "Random-forest consensus cleaning needs at least one label "
            "with two or more cells."
        )
    x_cv, y_cv = x.iloc[eligible], y.iloc[eligible]
    folds = min(cv_folds, int(y_cv.value_counts().min()))
    cv = RepeatedStratifiedKFold(
        n_splits=folds, n_repeats=repeats, random_state=seed
    )
    cleaner = RandomForestClassifier(
        n_estimators=100, max_features=mtry or "sqrt",
        n_jobs=num_threads, random_state=seed,
    )
    classes = np.unique(y_cv)
    probabilities = np.zeros((len(x), len(classes)))
    counts = np.zeros(len(x))
    match_counts = np.zeros(len(x))
    eligible_positions = np.flatnonzero(eligible)
    for train, test in cv.split(x_cv, y_cv):
        cleaner.fit(x_cv.iloc[train], y_cv.iloc[train])
        fold = cleaner.predict_proba(x_cv.iloc[test])
        positions = [np.flatnonzero(classes == value)[0] for value in cleaner.classes_]
        original_test = eligible_positions[test]
        probabilities[np.ix_(original_test, positions)] += fold
        counts[original_test] += 1
        fold_labels = cleaner.classes_[fold.argmax(axis=1)]
        match_counts[original_test] += (
            fold_labels == y_cv.iloc[test].to_numpy()
        )
    probabilities[eligible] /= counts[eligible, None]
    # R increments a match counter for every repeat's out-of-fold hard
    # prediction, then divides by repeats. Using only the winner of averaged
    # probabilities incorrectly collapses agreement to 0 or 1.
    rates = np.zeros(len(x))
    rates[eligible] = match_counts[eligible] / counts[eligible]
    retained = eligible & (rates >= agreement_thresh)
    cleaned = labels.copy()
    cleaned.loc[x.index[~retained]] = unknown_label
    spe.obs["cleaned_core_label"] = cleaned
    model = RandomForestClassifier(
        n_estimators=num_trees, max_features=mtry or "sqrt",
        n_jobs=num_threads, random_state=seed,
    ).fit(x.loc[retained], y.loc[retained])
    model.cytogater_features_ = features_used
    return {
        "spe": spe, "model": model,
        "core_prob_mat": pd.DataFrame(
            probabilities[retained], index=x.index[retained], columns=classes
        ),
        "agreement_rates": pd.Series(rates[retained], index=x.index[retained]),
        "features_used": features_used,
    }


def predict_unknown_with_randomforest(
    spe,
    model,
    assay_name="exprs",
    label_col="custom_label",
    out_col="soft_tree_label_filled",
    pred_col="rf_pred",
    unknown_label="Unknown",
    threshold=0.5,
    unassigned_label="Unassigned",
):
    frame = expression_frame(spe, assay_name)
    # AnnData commonly round-trips string observation columns as Categoricals.
    # Prediction may introduce a model class or ``Unassigned`` that is not an
    # existing category, so perform replacement on an object-backed Series.
    labels = spe.obs[label_col].astype(object).copy()
    replace = labels.isna() | (labels == unknown_label)
    if not replace.any():
        return {"spe": spe, "prob_mat": None}
    features = getattr(model, "cytogater_features_", list(frame.columns))
    values = model.predict_proba(frame.loc[replace, features])
    probabilities = pd.DataFrame(
        values, index=frame.index[replace], columns=model.classes_
    )
    confidence = probabilities.max(axis=1)
    predicted = probabilities.idxmax(axis=1).where(
        confidence >= threshold, unassigned_label
    )
    spe.obs[pred_col] = pd.Series(index=spe.obs_names, dtype=object)
    spe.obs.loc[replace, pred_col] = predicted
    labels.loc[replace] = predicted
    spe.obs[out_col] = labels
    spe.obs["rf_confidence"] = pd.Series(confidence, index=confidence.index)
    return {"spe": spe, "prob_mat": probabilities}


def rf_metric_table(fit):
    from ..metrics import class_metrics_from_fit
    class_table = class_metrics_from_fit(fit)
    return {"class": class_table}


def rf_metric_text(fit, digits=3):
    table = rf_metric_table(fit)["class"]
    return [
        f"{row['class']}: precision={row.precision:.{digits}f}, "
        f"recall={row.recall:.{digits}f}, f1={row.f1:.{digits}f}"
        for _, row in table.iterrows()
    ]
