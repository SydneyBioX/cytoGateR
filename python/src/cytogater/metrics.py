"""Label confidence and classification metrics."""

import numpy as np
import pandas as pd


def compute_custom_labels(
    prob_mat, base_thresh=0.4, quantile_cut=0.75, iqr_mult=1.5, flag_fn=None
):
    values = np.asarray(prob_mat, dtype=float)
    flags = []
    for row in values:
        finite = row[np.isfinite(row)]
        if flag_fn is not None:
            flags.append(bool(flag_fn(row)))
        elif not len(finite):
            flags.append(False)
        else:
            q = np.quantile(finite, quantile_cut)
            iqr = np.quantile(finite, 0.75) - np.quantile(finite, 0.25)
            flags.append(np.max(finite) > max(base_thresh, q + iqr_mult * iqr))
    return np.asarray(flags)


def assign_confident_labels(prob_mat, confidence_flags, unknown_label="Unknown"):
    values = np.asarray(prob_mat, dtype=float)
    columns = np.asarray(prob_mat.columns)
    safe = np.where(np.isfinite(values), values, -np.inf)
    labels = columns[safe.argmax(axis=1)].astype(object)
    labels[~np.isfinite(values).any(axis=1)] = unknown_label
    labels[~np.asarray(confidence_flags, dtype=bool)] = unknown_label
    return labels


def custom_labels(
    spe,
    prob_mat,
    colname="custom_label",
    base_thresh=0.4,
    quantile_cut=0.75,
    iqr_mult=1.5,
    unknown_label="Unknown",
    flag_fn=None,
):
    flags = compute_custom_labels(
        prob_mat, base_thresh, quantile_cut, iqr_mult, flag_fn
    )
    spe.obs[colname] = assign_confident_labels(prob_mat, flags, unknown_label)
    return spe


def class_metrics_from_fit(fit, test_truth=None, test_pred=None):
    if "cv_class_metrics" in fit:
        return fit["cv_class_metrics"]
    truth = test_truth if test_truth is not None else fit.get("test_truth")
    pred = test_pred if test_pred is not None else fit.get("test_pred")
    if truth is None or pred is None:
        raise ValueError("Provide test_truth/test_pred or metrics in fit.")
    truth, pred = np.asarray(truth).astype(str), np.asarray(pred).astype(str)
    rows = []
    for cls in sorted(set(truth) | set(pred)):
        tp = int(np.sum((truth == cls) & (pred == cls)))
        fp = int(np.sum((truth != cls) & (pred == cls)))
        fn = int(np.sum((truth == cls) & (pred != cls)))
        precision = tp / (tp + fp) if tp + fp else 0
        recall = tp / (tp + fn) if tp + fn else 0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0
        rows.append({
            "class": cls, "tp": tp, "fp": fp, "fn": fn,
            "support": int(np.sum(truth == cls)), "precision": precision,
            "recall": recall, "f1": f1,
        })
    return pd.DataFrame(rows)


def calculate_f1(spe, ref_col="ref_broad", pred_col="pred_broad"):
    table = class_metrics_from_fit(
        {}, spe.obs[ref_col].astype(str), spe.obs[pred_col].astype(str)
    )
    table = table[~table["class"].str.contains(
        "unassigned|undefined|unknown", case=False, regex=True
    )]
    return pd.DataFrame({
        "Category": table["class"],
        "Precision": table.precision.round(3),
        "Recall": table.recall.round(3),
        "F1_Score": table.f1.round(3),
        "Cell_Count": table.support,
    }).reset_index(drop=True)


def label_confusion_matrix(spe, label_col1, label_col2, drop_na=True):
    frame = spe.obs[[label_col1, label_col2]]
    if drop_na:
        frame = frame.dropna()
    levels = sorted(set(frame[label_col1].astype(str)) | set(frame[label_col2].astype(str)))
    return pd.crosstab(
        pd.Categorical(frame[label_col1], levels),
        pd.Categorical(frame[label_col2], levels),
        dropna=False,
    )


def label_agreement_rates(spe, label_col1, label_col2, drop_na=True):
    frame = spe.obs[[label_col1, label_col2]]
    if drop_na:
        frame = frame.dropna()
    a, b = frame[label_col1].astype(str), frame[label_col2].astype(str)
    rows = []
    for label in sorted(set(a) | set(b)):
        match = int(((a == label) & (b == label)).sum())
        union = int(((a == label) | (b == label)).sum())
        rows.append({
            "label": label, "match_n": match, "union_n": union,
            "agreement": match / union if union else np.nan,
        })
    return pd.DataFrame(rows)


def prob_quantile_cutoff(x, prob=0.98):
    return float(np.nanquantile(x, prob))


def prob_mad_cutoff(x, mad_mult=3):
    x = np.asarray(x, dtype=float)
    median = np.nanmedian(x)
    # R uses stats::mad(..., constant = 1) in prob_mad_cutoff().
    mad = np.nanmedian(np.abs(x - median))
    return float(median + mad_mult * mad)


def probability_label_matrix(prob_mat, cutoff_fn=prob_quantile_cutoff):
    frame = pd.DataFrame(prob_mat)
    cutoffs = frame.apply(lambda column: cutoff_fn(column))
    cutoffs[~np.isfinite(cutoffs)] = np.nan
    label_matrix = frame.ge(cutoffs, axis="columns").astype("boolean")
    missing = frame.isna()
    missing.loc[:, cutoffs.isna()] = True
    return label_matrix.mask(missing)


def apply_cutoff_labels(
    res,
    cutoff_fn=prob_quantile_cutoff,
    label_col="cutoff_label",
    unknown_label="Unknown",
):
    labels = probability_label_matrix(res["prob_mat"], cutoff_fn)
    positive_counts = labels.sum(axis=1, skipna=True)
    names = pd.Series(unknown_label, index=labels.index, dtype=object)
    single = positive_counts == 1
    if single.any():
        names.loc[single] = labels.loc[single].fillna(False).idxmax(axis=1)
    res["spe"].obs[label_col] = names.reindex(res["spe"].obs_names).to_numpy()
    res["label_mat"] = labels
    return res
