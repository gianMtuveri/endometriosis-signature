from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from endometriosis_signature.dataset import get_feature_columns


def make_l1_classifier(C: float) -> Pipeline:
    return Pipeline([
        ("scaler", StandardScaler()),
        ("model", LogisticRegression(
            penalty="l1",
            solver="liblinear",
            C=C,
            max_iter=2000,
        )),
    ])


def main() -> None:
    processed_dir = Path("data/processed")
    processed_dir.mkdir(parents=True, exist_ok=True)

    # =========================
    # LOAD DATA
    # =========================
    data = pd.read_parquet(processed_dir / "endometriosis_clean_filtered.parquet")
    feature_cols = get_feature_columns(data)

    X = data[feature_cols]
    y = data["target"]

    # =========================
    # LOAD RANKED PROBES
    # =========================
    ranked = pd.read_csv(processed_dir / "ranked_signature_probes.csv")

    ranked = ranked[
        ranked["Gene Symbol"].notna()
        & (ranked["Gene Symbol"].astype(str).str.strip() != "")
    ].copy()

    # keep the current ranked order from your main pipeline
    top_probes = ranked["probe_id"].tolist()

    # use the compact ranked subset, not the whole transcriptome
    # this makes the L1 test focused and interpretable
    X_ranked = X[top_probes]

    print("Input matrix for L1:", X_ranked.shape)

    # =========================
    # SWEEP C VALUES
    # =========================
    c_values = [0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10]
    rows = []

    for C in c_values:
        clf = make_l1_classifier(C)

        scores = cross_val_score(
            clf,
            X_ranked,
            y,
            cv=5,
            scoring="roc_auc"
        )

        # fit on full data to inspect sparsity
        clf.fit(X_ranked, y)
        coef = clf.named_steps["model"].coef_.flatten()

        nonzero_mask = coef != 0
        n_nonzero = int(nonzero_mask.sum())

        selected_probes = np.array(top_probes)[nonzero_mask].tolist()

        # map to genes when possible
        selected_df = ranked[ranked["probe_id"].isin(selected_probes)][
            ["probe_id", "Gene Symbol", "Gene Title"]
        ].drop_duplicates()

        selected_genes = selected_df["Gene Symbol"].tolist()

        rows.append({
            "C": C,
            "mean_auc": float(scores.mean()),
            "std_auc": float(scores.std()),
            "n_nonzero": n_nonzero,
            "selected_probes": ",".join(selected_probes),
            "selected_genes": ",".join(selected_genes),
        })

        print(
            f"C={C:<4} | AUC={scores.mean():.3f} ± {scores.std():.3f} | "
            f"non-zero={n_nonzero}"
        )

    results = pd.DataFrame(rows)
    results.to_csv(processed_dir / "signature_l1_sweep.csv", index=False)

    print("\nSaved results to data/processed/signature_l1_sweep.csv")

    # =========================
    # PRINT BEST ROWS
    # =========================
    best = results.sort_values(
        ["mean_auc", "n_nonzero"],
        ascending=[False, True]
    ).head(5)

    print("\nTop L1 configurations:")
    print(best[["C", "mean_auc", "std_auc", "n_nonzero", "selected_genes"]])


if __name__ == "__main__":
    main()
