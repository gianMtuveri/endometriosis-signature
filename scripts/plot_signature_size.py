from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from endometriosis_signature.dataset import get_feature_columns
from endometriosis_signature.modeling import make_classifier


processed_dir = Path("data/processed")

data = pd.read_parquet(processed_dir / "endometriosis_clean_filtered.parquet")
feature_cols = get_feature_columns(data)

X = data[feature_cols]
y = data["target"]

ranked = pd.read_csv(processed_dir / "ranked_signature_probes.csv")

ranked = ranked[
    ranked["Gene Symbol"].notna()
    & (ranked["Gene Symbol"].astype(str).str.strip() != "")
].copy()

top_probes = ranked["probe_id"].tolist()

print(f"Number of ranked probes available: {len(top_probes)}")
print("Top probes:")
print(ranked[["probe_id", "count", "coef", "Gene Symbol"]].head(10))

sizes = list(range(2, len(top_probes) + 6))
means = []
stds = []

for k in sizes:
    selected = top_probes[:k]
    X_k = X[selected]

    clf = make_classifier()

    from sklearn.model_selection import cross_val_score
    scores = cross_val_score(
        clf,
        X_k,
        y,
        cv=5,
        scoring="roc_auc"
    )

    means.append(scores.mean())
    stds.append(scores.std())

    print(f"{k} genes -> AUC: {scores.mean():.3f} ± {scores.std():.3f}")

plt.figure(figsize=(7, 5))
plt.plot(sizes, means, marker="o", linewidth=2)
plt.fill_between(
    sizes,
    np.array(means) - np.array(stds),
    np.array(means) + np.array(stds),
    alpha=0.2
)

best_idx = int(np.argmax(means))
best_k = sizes[best_idx]
plt.axvline(best_k, linestyle="--")

plt.xlabel("Number of genes in signature")
plt.ylabel("ROC-AUC")
plt.title("Signature Size vs Performance")
plt.tight_layout()

out = "results/figures/signature_size_curve.png"
plt.savefig(out, dpi=300)
plt.show()

print(f"Saved figure to: {out}")
print(f"Best k by mean AUC: {best_k}")