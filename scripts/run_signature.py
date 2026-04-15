from pathlib import Path

import pandas as pd

from endometriosis_signature.annotation import annotate_probes, load_gpl570_annotation
from endometriosis_signature.dataset import get_feature_columns
from endometriosis_signature.modeling import (
    compute_stable_features,
    evaluate_auc_cv,
    fit_and_rank_features,
    make_classifier,
    signature_size_sweep,
)


def main() -> None:
    processed_dir = Path("data/processed")
    raw_dir = Path("data/raw")

    data = pd.read_parquet(processed_dir / "endometriosis_clean_filtered.parquet")
    feature_cols = get_feature_columns(data)

    X = data[feature_cols]
    y = data["target"]

    clf = make_classifier()

    baseline = evaluate_auc_cv(clf, X, y)
    print(f"Baseline AUC: {baseline['mean_auc']:.3f} +- {baseline['std_auc']:.3f}")

    stable_df = compute_stable_features(clf, X, y, feature_cols, n_splits=5, top_n=20)
    stable_df.to_csv(processed_dir / "stable_probe_counts.csv", index=False)

    gpl_small = load_gpl570_annotation(raw_dir / "GPL570-tbl-1.txt")
    annotated = annotate_probes(stable_df, gpl_small)

    stable_mapped = annotated[
        (annotated["count"] >= 3)
        & (annotated["Gene Symbol"].notna())
        & (annotated["Gene Symbol"].astype(str).str.strip() != "")
    ].copy()
    stable_mapped.to_csv(processed_dir / "stable_mapped_probes.csv", index=False)

    ranked = fit_and_rank_features(clf, X, y, feature_cols, stable_mapped)
    ranked.to_csv(processed_dir / "ranked_signature_probes.csv", index=False)

    print("\nTop ranked probes:")
    print(ranked[["probe_id", "count", "coef", "Gene Symbol", "Gene Title"]].head(10))

    sweep_df = signature_size_sweep(clf, data, y, ranked, k_min=2, k_max=14)
    sweep_df.to_csv(processed_dir / "signature_size_sweep.csv", index=False)

    print("\nSignature size sweep:")
    for _, row in sweep_df.iterrows():
        print(f"{int(row['k']):2d} genes -> AUC: {row['mean_auc']:.3f} ± {row['std_auc']:.3f}")

    best_row = sweep_df.sort_values(["mean_auc", "std_auc"], ascending=[False, True]).iloc[0]
    best_k = int(best_row["k"])
    best_signature = ranked.head(best_k)["probe_id"].tolist()

    print(f"\nBest signature size: {best_k}")
    print("Best signature probes:", best_signature)


if __name__ == "__main__":
    main()