import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import StratifiedKFold

# =========================
# SIGNATURE DEFINITION
# =========================

gene_map = {
    "CTU2": "1561501_s_at",
    "ZNF24": "242210_at",
    "NT5DC3": "234737_at",
    "HMGN3-AS1": "1559404_a_at",
    "ZNF568": "1560779_a_at",
    "C11orf54": "1559623_at",
}

genes = list(gene_map.keys())
probes = list(gene_map.values())

# =========================
# LOAD TRAIN DATA
# =========================

data = pd.read_parquet("data/processed/endometriosis_clean_filtered.parquet")

X_train = data[probes].copy()
X_train.columns = genes
y_train = data["target"]

# =========================
# LOAD EXTERNAL DATA
# =========================

df_new = pd.read_csv(
    "data/raw/GSE135485_Endometriosis_raw_counts.csv",
    index_col=0
)

# ---- extract labels ----
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

y_test = pd.Series(
    [
        1 if "patient with endometriosis" in x.lower() else 0
        for x in status_row
    ],
    index=df_new.columns,
)

# ---- match genes ----
available = [g for g in genes if g in df_new.index]

X_test = np.log1p(df_new.loc[available].T)

# =========================
# MODEL
# =========================

clf = Pipeline([
    ("scaler", StandardScaler()),
    ("model", LogisticRegression(max_iter=2000))
])

# =========================
# TRAIN ROC (CV)
# =========================

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

y_true_all = []
y_score_all = []

for train_idx, test_idx in skf.split(X_train[available], y_train):
    X_tr = X_train.iloc[train_idx][available]
    y_tr = y_train.iloc[train_idx]

    X_val = X_train.iloc[test_idx][available]
    y_val = y_train.iloc[test_idx]

    clf.fit(X_tr, y_tr)
    proba = clf.predict_proba(X_val)[:, 1]

    y_true_all.extend(y_val)
    y_score_all.extend(proba)

fpr_train, tpr_train, _ = roc_curve(y_true_all, y_score_all)
auc_train = roc_auc_score(y_true_all, y_score_all)

# =========================
# TEST ROC
# =========================

clf.fit(X_train[available], y_train)
proba_test = clf.predict_proba(X_test[available])[:, 1]

fpr_test, tpr_test, _ = roc_curve(y_test, proba_test)
auc_test = roc_auc_score(y_test, proba_test)

# =========================
# PLOT
# =========================

plt.figure()

plt.plot(fpr_train, tpr_train, label=f"Train (CV) AUC = {auc_train:.2f}")
plt.plot(fpr_test, tpr_test, label=f"External AUC = {auc_test:.2f}")

plt.plot([0, 1], [0, 1], linestyle="--")

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve - Endometriosis Signature")
plt.legend()

# save figure
plt.savefig("results/figures/roc_curve.png", dpi=300)
plt.show()
