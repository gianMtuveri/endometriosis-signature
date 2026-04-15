from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd
from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def make_classifier() -> Pipeline:
    return Pipeline(
        [
            ("var", VarianceThreshold(threshold=0.01)),
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(max_iter=2000)),
        ]
    )


def evaluate_auc_cv(clf: Pipeline, X: pd.DataFrame, y: pd.Series, cv: int = 5) -> dict:
    scores = cross_val_score(clf, X, y, cv=cv, scoring="roc_auc")
    return {
        "scores": scores,
        "mean_auc": float(scores.mean()),
        "std_auc": float(scores.std()),
    }


def compute_stable_features(
    clf: Pipeline,
    X: pd.DataFrame,
    y: pd.Series,
    feature_cols: list[str],
    n_splits: int = 5,
    top_n: int = 20,
) -> pd.DataFrame:
    feature_counts = Counter()
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    for train_idx, _ in skf.split(X, y):
        X_train = X.iloc[train_idx]
        y_train = y.iloc[train_idx]

        clf.fit(X_train, y_train)

        var = clf.named_steps["var"]
        model = clf.named_steps["model"]

        mask = var.get_support()
        selected_features = np.array(feature_cols)[mask]

        coef = model.coef_.flatten()
        top_idx = np.argsort(np.abs(coef))[-top_n:]
        top_features = selected_features[top_idx]

        feature_counts.update(top_features)

    stable_df = pd.DataFrame(
        sorted(feature_counts.items(), key=lambda x: (-x[1], x[0])),
        columns=["probe_id", "count"],
    )
    return stable_df


def fit_and_rank_features(
    clf: Pipeline,
    X: pd.DataFrame,
    y: pd.Series,
    feature_cols: list[str],
    stable_mapped: pd.DataFrame,
) -> pd.DataFrame:
    clf.fit(X, y)

    var = clf.named_steps["var"]
    model = clf.named_steps["model"]

    mask = var.get_support()
    selected_features = np.array(feature_cols)[mask]
    coef = model.coef_.flatten()

    coef_df = pd.DataFrame({"probe_id": selected_features, "coef": coef})
    coef_df["abs_coef"] = coef_df["coef"].abs()

    ranked = stable_mapped.merge(coef_df, on="probe_id", how="left")
    ranked = ranked.sort_values(["count", "abs_coef"], ascending=[False, False])
    return ranked


def signature_size_sweep(
    clf: Pipeline,
    data: pd.DataFrame,
    y: pd.Series,
    ranked: pd.DataFrame,
    k_min: int = 2,
    k_max: int = 14,
) -> pd.DataFrame:
    results = []

    for k in range(k_min, k_max + 1):
        sig = ranked.head(k)["probe_id"].tolist()
        X_sig = data[sig]

        scores_k = cross_val_score(clf, X_sig, y, cv=5, scoring="roc_auc")
        results.append(
            {
                "k": k,
                "mean_auc": float(scores_k.mean()),
                "std_auc": float(scores_k.std()),
            }
        )

    return pd.DataFrame(results)