import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score

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

# =========================
# 1. LOAD GSE25628 MATRIX
# =========================

df = pd.read_csv(
    "data/raw/GSE25628_series_matrix.txt",
    sep="\t",
    comment="!",
    index_col=0,
)

print("Raw GSE25628 shape:", df.shape)

# =========================
# 2. EXTRACT LABELS
#    Keep only:
#    1 -> Pathological (ectopic)
#    0 -> Normal (control)
#    Exclude:
#    Normal (eutopic)
# =========================

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

y_test = pd.Series(labels, index=df.columns, name="target")
mask = y_test.notna()

# Keep only selected samples
df = df.loc[:, mask]
y_test = y_test[mask].astype(int)

print("\nFiltered GSE25628 shape:", df.shape)
print("Filtered label counts:")
print(y_test.value_counts())

# =========================
# 3. LOAD GPL571 ANNOTATION
# =========================

gpl = pd.read_csv(
    "data/raw/GPL571-tbl-1.txt",
    sep="\t",
    comment="!",
    header=None,
    names=GPL571_COLUMNS,
    low_memory=False,
)

gpl_small = gpl[["ID", "Gene Symbol", "Gene Title"]].copy()
gpl_small = gpl_small.rename(columns={"ID": "probe_id"})

# =========================
# 4. MAP GSE25628 PROBES -> GENE SYMBOLS
# =========================

df_reset = df.reset_index().rename(columns={"ID_REF": "probe_id"})
df_annot = df_reset.merge(gpl_small, on="probe_id", how="left")

df_annot = df_annot[df_annot["Gene Symbol"].notna()].copy()
df_annot["Gene Symbol"] = df_annot["Gene Symbol"].astype(str).str.strip()

# If multiple probes map to the same gene, average them
sample_cols = list(df.columns)
df_gene = df_annot.groupby("Gene Symbol")[sample_cols].mean()

print("\nGene-level GSE25628 shape:", df_gene.shape)

# =========================
# 5. DEFINE TRANSFERABLE SIGNATURE
# =========================

signature_genes = [
    "CTU2",
    "ZNF24",
    "NT5DC3",
    "HMGN3-AS1",
    "ZNF568",
    "C11orf54",
]

available = [g for g in signature_genes if g in df_gene.index]
print("\nAvailable signature genes in GSE25628:", available)

if len(available) == 0:
    raise ValueError("No signature genes available in GSE25628 after annotation.")

X_test = df_gene.loc[available, y_test.index].T
print("External X_test shape:", X_test.shape)

# =========================
# 6. LOAD TRAINING DATASET
#    GSE51981 clean subset
# =========================

train_data = pd.read_parquet("data/processed/endometriosis_clean_filtered.parquet")

gene_map = {
    "CTU2": "1561501_s_at",
    "ZNF24": "242210_at",
    "NT5DC3": "234737_at",
    "HMGN3-AS1": "1559404_a_at",
    "ZNF568": "1560779_a_at",
    "C11orf54": "1559623_at",
}

X_train = train_data[[gene_map[g] for g in signature_genes]].copy()
X_train.columns = signature_genes
y_train = train_data["target"]

X_train = X_train[available]

print("Training X_train shape:", X_train.shape)

# =========================
# 7. TRAIN ON GSE51981
# =========================

clf = Pipeline([
    ("scaler", StandardScaler()),
    ("model", LogisticRegression(max_iter=2000))
])

clf.fit(X_train, y_train)

# =========================
# 8. TEST ON GSE25628
# =========================

proba = clf.predict_proba(X_test[available])[:, 1]

auc = roc_auc_score(y_test, proba)
pr = average_precision_score(y_test, proba)

print("\n=== EXTERNAL VALIDATION: GSE25628 ===")
print("ROC-AUC:", auc)
print("PR-AUC:", pr)

# =========================
# 9. SAVE RESULTS
# =========================

results = pd.DataFrame([
    {
        "dataset": "GSE25628",
        "platform": "GPL571",
        "n_genes_used": len(available),
        "genes_used": ",".join(available),
        "n_samples": len(y_test),
        "n_positive": int((y_test == 1).sum()),
        "n_negative": int((y_test == 0).sum()),
        "roc_auc": auc,
        "pr_auc": pr,
        "notes": "Pathological (ectopic) vs Normal (control); eutopic excluded",
    }
])

out_path = "data/processed/external_validation_gse25628.csv"
results.to_csv(out_path, index=False)

print(f"\nSaved results to: {out_path}")