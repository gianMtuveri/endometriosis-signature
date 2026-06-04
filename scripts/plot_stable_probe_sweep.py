from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import cross_val_score

from endometriosis_signature.dataset import get_feature_columns
from endometriosis_signature.modeling import make_classifier


processed_dir = Path("data/processed")
figures_dir = Path("results/figures")
figures_dir.mkdir(parents=True, exist_ok=True)

data = pd.read_parquet(processed_dir / "endometriosis_clean_filtered.parquet")
feature_cols = get_feature_columns(data)

X = data[feature_cols]
y = data["target"]

stable = pd.read_csv(processed_dir / "stable_probe_counts.csv")

# Keep only probes that are actually present in the expression matrix
stable = stable[stable["probe_id"].isin(feature_cols)].copy()

# Rank by stability count
stable = stable.sort_values("count", ascending=False)

top_probes = stable["probe_id"].tolist()

print(f"Available stable probes: {len(top_probes)}")
print(stable.head(20))

sizes = list(range(2, len(top_probes) + 1))

means = []
stds = []

for k in sizes:
    selected = top_probes[:k]
    X_k = X[selected]

    clf = make_classifier()

    scores = cross_val_score(
        clf,
        X_k,
        y,
        cv=5,
        scoring="roc_auc"
    )

    means.append(scores.mean())
    stds.append(scores.std())

    print(f"{k} probes -> AUC: {scores.mean():.3f} ± {scores.std():.3f}")

sweep = pd.DataFrame({
    "k": sizes,
    "mean_auc": means,
    "std_auc": stds,
})

sweep.to_csv(processed_dir / "stable_probe_size_sweep.csv", index=False)

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
plt.xlabel("Number of stable probes included")
plt.ylabel("ROC-AUC")
plt.title("Exploratory Stable-Probe Sweep")
plt.tight_layout()

out = figures_dir / "stable_probe_size_curve.png"
plt.savefig(out, dpi=300)
plt.show()

print(f"Saved figure to: {out}")
print(f"Best k by mean AUC: {best_k}")