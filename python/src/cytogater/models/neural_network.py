"""Optional PyTorch classifier with the R package's public workflow."""

import numpy as np
import pandas as pd

from ..preprocessing import expression_frame


def train_custom_dl(
    spe,
    label_col="cutoff_label",
    assay_name="exprs",
    unknown_label="Unknown",
    features="all",
    hidden_dims=(64, 32),
    dropout=0.3,
    epochs=50,
    lr=1e-3,
    weight_decay=1e-4,
    batch_size=256,
    cv_folds=5,
    repeats=3,
    agreement_thresh=0.8,
    val_frac=0.15,
    patience=10,
    verbose_training=False,
    class_balanced_loss=True,
    seed=None,
    **kwargs,
):
    """Train a small probability network; PyTorch is imported only when used."""
    import torch
    from torch import nn

    if seed is not None:
        torch.manual_seed(seed)
    frame = expression_frame(spe, assay_name)
    features_used = (
        list(frame.columns)
        if isinstance(features, str) and features == "all"
        else [x for x in features if x in frame]
    )
    labels = spe.obs[label_col]
    keep = labels.notna() & (labels != unknown_label)
    x = frame.loc[keep, features_used].to_numpy(np.float32)
    classes, y = np.unique(labels[keep].astype(str), return_inverse=True)
    mean, scale = x.mean(0), x.std(0)
    scale[scale == 0] = 1
    x = (x - mean) / scale
    layers, width = [], x.shape[1]
    for hidden in hidden_dims:
        layers += [nn.Linear(width, hidden), nn.ReLU(), nn.Dropout(dropout)]
        width = hidden
    layers.append(nn.Linear(width, len(classes)))
    del cv_folds, repeats, agreement_thresh, val_frac, patience
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = nn.Sequential(*layers).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    if class_balanced_loss:
        counts = np.bincount(y)
        weights = len(y) / (len(counts) * counts)
        loss_fn = nn.CrossEntropyLoss(
            weight=torch.tensor(weights, dtype=torch.float32, device=device)
        )
    else:
        loss_fn = nn.CrossEntropyLoss()
    dataset = torch.utils.data.TensorDataset(
        torch.tensor(x), torch.tensor(y, dtype=torch.long)
    )
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True)
    model.train()
    for _ in range(epochs):
        for xb, yb in loader:
            optimizer.zero_grad()
            loss = loss_fn(model(xb.to(device)), yb.to(device))
            loss.backward()
            optimizer.step()
        if verbose_training:
            print(f"epoch {_ + 1}/{epochs}: loss={loss.item():.4f}")
    fitted = {
        "net": model, "class_levels": classes, "feature_names": features_used,
        "center": mean, "scale": scale, "device": device,
        "training_options": kwargs,
    }
    return {"spe": spe, "model": fitted, "features_used": features_used}


def predict_unknown_with_dl(
    spe,
    model,
    assay_name="exprs",
    label_col="custom_label",
    out_col="soft_tree_label_filled",
    pred_col="dl_pred",
    unknown_label="Unknown",
    threshold=0.5,
    unassigned_label="Unassigned",
):
    import torch

    frame = expression_frame(spe, assay_name)
    labels = spe.obs[label_col].copy()
    replace = labels.isna() | (labels == unknown_label)
    if not replace.any():
        return {"spe": spe, "prob_mat": None}
    x = frame.loc[replace, model["feature_names"]].to_numpy(np.float32)
    x = (x - model["center"]) / model["scale"]
    model["net"].eval()
    with torch.no_grad():
        probability = torch.softmax(
            model["net"](torch.tensor(x).to(model["device"])), dim=1
        ).cpu().numpy()
    probs = pd.DataFrame(
        probability, index=frame.index[replace], columns=model["class_levels"]
    )
    confidence = probs.max(axis=1)
    predicted = probs.idxmax(axis=1).where(confidence >= threshold, unassigned_label)
    spe.obs[pred_col] = pd.Series(index=spe.obs_names, dtype=object)
    spe.obs.loc[replace, pred_col] = predicted
    labels.loc[replace] = predicted
    spe.obs[out_col] = labels
    spe.obs["dl_confidence"] = pd.Series(confidence, index=confidence.index)
    return {"spe": spe, "prob_mat": probs}
