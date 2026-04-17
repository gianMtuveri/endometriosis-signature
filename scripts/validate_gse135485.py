import pandas as pd
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score

DATASET_NAME = "gse135485"
SIGNATURE_NAME = "sig3"

# map genes -> probes (from your signature)
'''gene_map = { # genes map for 7 gene signature
    "CTU2": "1561501_s_at",
    "ZNF24": "242210_at",
    "NT5DC3": "234737_at",
    "HMGN3-AS1": "1559404_a_at",
    "ZNF568": "1560779_a_at",
    "C11orf54": "1559623_at",
}'''
gene_map = { # genes map for 3 gene signature
    "ZNF24": "242210_at",
    "HMGN3-AS1": "1559404_a_at",
    "ZNF568": "1560779_a_at",
}

# =========================
# 1. LOAD DATASET 1 (TRAIN)
# =========================


data = pd.read_parquet("data/processed/endometriosis_clean_filtered.parquet")


genes = list(gene_map.keys())
probes = list(gene_map.values())

# training data
X_train = data[probes].copy()
X_train.columns = genes
y_train = data["target"]

print("Train shape:", X_train.shape)

# =========================
# 2. LOAD DATASET 2 (TEST)
# =========================

df_new = pd.read_csv(
    "data/raw/GSE135485_Endometriosis_raw_counts.csv",
    index_col=0
)

print("New dataset shape:", df_new.shape)

# =========================
# 3. EXTRACT LABELS
# =========================

meta_lines = []
with open("data/raw/GSE135485_series_matrix.txt", "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("!Sample_"):
            meta_lines.append(line.rstrip("\n").split("\t"))

status_row = None
for row in meta_lines:
    if row[0] == "!Sample_characteristics_ch1":
        values = [x.strip('"') for x in row[1:]]
        if any("subject status:" in x.lower() for x in values):
            status_row = values
            break

if status_row is None:
    raise ValueError("Could not find subject status")

y_new = pd.Series(
    [
        1 if "patient with endometriosis" in x.lower() else 0
        for x in status_row
    ],
    index=df_new.columns,
    name="target",
)

print("\nExternal label distribution:")
print(y_new.value_counts())

# =========================
# 4. MATCH GENES
# =========================

available = [g for g in genes if g in df_new.index]

print("\nAvailable genes in dataset 2:", available)

# =========================
# 5. BUILD TEST MATRIX
# =========================

X_new = np.log1p(df_new.loc[available].T)

print("External X shape:", X_new.shape)

# =========================
# 6. TRAIN MODEL (ON DATASET 1)
# =========================

clf = Pipeline([
    ("scaler", StandardScaler()),
    ("model", LogisticRegression(max_iter=2000))
])

X_train_common = X_train[available]
clf.fit(X_train_common, y_train)

# =========================
# 7. TEST ON DATASET 2
# =========================

proba = clf.predict_proba(X_new[available])[:, 1]

# =========================
# 8. METRICS
# =========================

auc = roc_auc_score(y_new, proba)
pr = average_precision_score(y_new, proba)

print("\n=== EXTERNAL VALIDATION ===")
print("ROC-AUC:", auc)
print("PR-AUC:", pr)

# =========================
# 9. RESULTS
# =========================

results = pd.DataFrame([
    {
        "n_genes": len(available),
        "genes": ",".join(available),
        "roc_auc": auc,
        "pr_auc": pr,
        "n_samples": len(y_new),
        "n_positive": int((y_new == 1).sum()),
        "n_negative": int((y_new == 0).sum()),
        "notes": "Cross-platform (microarray→RNA-seq); heavy class imbalance"
    }
])

out_path = f"data/processed/validation_{DATASET_NAME}_{SIGNATURE_NAME}.csv"
results.to_csv(out_path, index=False)

print(f"\nSaved external results to: {out_path}")
