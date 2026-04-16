import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import StratifiedKFold


# =========================
# CONFIG
# =========================

gene_map = {
    "CTU2": "1561501_s_at",
    "ZNF24": "242210_at",
    "NT5DC3": "234737_at",
    "HMGN3-AS1": "1559404_a_at",
    "ZNF568": "1560779_a_at",
    "C11orf54": "1559623_at",
}
signature_genes = list(gene_map.keys())


def make_model() -> Pipeline:
    return Pipeline([
        ("scaler", StandardScaler()),
        ("model", LogisticRegression(max_iter=2000)),
    ])


# =========================
# 1. TRAIN DATA
# =========================

train_data = pd.read_parquet("data/processed/endometriosis_clean_filtered.parquet")

X_train_full = train_data[[gene_map[g] for g in signature_genes]].copy()
X_train_full.columns = signature_genes
y_train = train_data["target"]

clf = make_model()

# internal CV ROC
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

y_true_cv = []
y_score_cv = []

for train_idx, test_idx in skf.split(X_train_full, y_train):
    X_tr = X_train_full.iloc[train_idx]
    y_tr = y_train.iloc[train_idx]

    X_val = X_train_full.iloc[test_idx]
    y_val = y_train.iloc[test_idx]

    clf.fit(X_tr, y_tr)
    proba = clf.predict_proba(X_val)[:, 1]

    y_true_cv.extend(y_val)
    y_score_cv.extend(proba)

fpr_train, tpr_train, _ = roc_curve(y_true_cv, y_score_cv)
auc_train = roc_auc_score(y_true_cv, y_score_cv)


# =========================
# 2. EXTERNAL RNA-SEQ (GSE135485)
# =========================

df_rna = pd.read_csv(
    "data/raw/GSE135485_Endometriosis_raw_counts.csv",
    index_col=0,
)

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
    raise ValueError("Could not find subject status row in GSE135485 metadata.")

y_rna = pd.Series(
    [1 if "patient with endometriosis" in x.lower() else 0 for x in status_row],
    index=df_rna.columns,
    name="target",
)

available_rna = [g for g in signature_genes if g in df_rna.index]

X_rna = np.log1p(df_rna.loc[available_rna].T)

X_train_rna = X_train_full[available_rna]
clf_rna = make_model()
clf_rna.fit(X_train_rna, y_train)

proba_rna = clf_rna.predict_proba(X_rna[available_rna])[:, 1]

fpr_rna, tpr_rna, _ = roc_curve(y_rna, proba_rna)
auc_rna = roc_auc_score(y_rna, proba_rna)


# =========================
# 3. EXTERNAL MICROARRAY (GSE25628)
# =========================

GPL571_COLUMNS = [
    "ID",
    "GB_ACC",
    "SPOT_ID",
    "Species Scientific Name",
    "Annotation Date",
    "Sequence Type",
    "Sequence Source",
    "Target Description",
    "Representative Public ID",
    "Gene Title",
    "Gene Symbol",
    "ENTREZ_GENE_ID",
    "RefSeq Transcript ID",
    "Gene Ontology Biological Process",
    "Gene Ontology Cellular Component",
    "Gene Ontology Molecular Function",
]

df_ma = pd.read_csv(
    "data/raw/GSE25628_series_matrix.txt",
    sep="\t",
    comment="!",
    index_col=0,
)

meta_lines = []
with open("data/raw/GSE25628_series_matrix.txt", "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("!Sample_"):
            meta_lines.append(line.rstrip("\n").split("\t"))

labels_row = None
for row in meta_lines:
    if row[0] == "!Sample_characteristics_ch1":
        values = [x.strip('"').lower() for x in row[1:]]
        if any("disease state" in x for x in values):
            labels_row = values
            break

if labels_row is None:
    raise ValueError("Could not find disease state row in GSE25628 metadata.")

labels = []
for x in labels_row:
    if "pathological (ectopic)" in x:
        labels.append(1)
    elif "normal (control)" in x:
        labels.append(0)
    else:
        labels.append(None)  # exclude eutopic

y_ma = pd.Series(labels, index=df_ma.columns, name="target")
mask = y_ma.notna()

df_ma = df_ma.loc[:, mask]
y_ma = y_ma[mask].astype(int)

gpl571 = pd.read_csv(
    "data/raw/GPL571-tbl-1.txt",
    sep="\t",
    comment="!",
    header=None,
    names=GPL571_COLUMNS,
    low_memory=False,
)

gpl571_small = gpl571[["ID", "Gene Symbol"]].copy()
gpl571_small = gpl571_small.rename(columns={"ID": "probe_id"})

df_ma_reset = df_ma.reset_index().rename(columns={"ID_REF": "probe_id"})
df_ma_annot = df_ma_reset.merge(gpl571_small, on="probe_id", how="left")
df_ma_annot = df_ma_annot[df_ma_annot["Gene Symbol"].notna()].copy()
df_ma_annot["Gene Symbol"] = df_ma_annot["Gene Symbol"].astype(str).str.strip()

sample_cols = list(df_ma.columns)
df_ma_gene = df_ma_annot.groupby("Gene Symbol")[sample_cols].mean()

available_ma = [g for g in signature_genes if g in df_ma_gene.index]

X_ma = df_ma_gene.loc[available_ma, y_ma.index].T

X_train_ma = X_train_full[available_ma]
clf_ma = make_model()
clf_ma.fit(X_train_ma, y_train)

proba_ma = clf_ma.predict_proba(X_ma[available_ma])[:, 1]

fpr_ma, tpr_ma, _ = roc_curve(y_ma, proba_ma)
auc_ma = roc_auc_score(y_ma, proba_ma)


# =========================
# 4. PLOT
# =========================

plt.figure(figsize=(7, 6))

plt.plot(fpr_train, tpr_train, linewidth=2, label=f"Internal CV (AUC = {auc_train:.2f})")
plt.plot(fpr_rna, tpr_rna, linewidth=2, label=f"GSE135485 RNA-seq (AUC = {auc_rna:.2f})")
plt.plot(fpr_ma, tpr_ma, linewidth=2, label=f"GSE25628 microarray (AUC = {auc_ma:.2f})")

plt.plot([0, 1], [0, 1], linestyle="--", linewidth=1)

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curves for Endometriosis Signature")
plt.legend(loc="lower right")
plt.tight_layout()

out_path = "results/figures/roc_curve_all.png"
plt.savefig(out_path, dpi=300)
plt.show()

print(f"Saved figure to: {out_path}")
print(f"Internal CV AUC: {auc_train:.3f}")
print(f"GSE135485 AUC: {auc_rna:.3f}")
print(f"GSE25628 AUC: {auc_ma:.3f}")
