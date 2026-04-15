import pandas as pd
import numpy as np
import json
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from collections import Counter

gpl_cols = [
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


data = pd.read_parquet("data/processed/endometriosis_clean_filtered.parquet")

meta_cols = ["target", "source", "tissue", "severity"]
feature_cols = [col for col in data.columns if col not in meta_cols]

X = data[feature_cols]
y = data["target"]

clf = Pipeline([
    ("var", VarianceThreshold(threshold=0.01)),
    ("scaler", StandardScaler()),
    ("model", LogisticRegression(max_iter=2000))
])

scores = cross_val_score(clf, X, y, cv=5, scoring="roc_auc")

#print("Baseline AUC:", scores.mean(), "+-", scores.std())

feature_counts = Counter()
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

for train_idx, test_idx in skf.split(X, y):
    X_train = X.iloc[train_idx]
    y_train = y.iloc[train_idx]

    clf.fit(X_train, y_train)

    var = clf.named_steps["var"]
    model = clf.named_steps["model"]

    mask = var.get_support()
    selected_features = np.array(feature_cols)[mask]

    coef = model.coef_.flatten()
    top_idx = np.argsort(np.abs(coef))[-20:]
    top_features = selected_features[top_idx]

    feature_counts.update(top_features)


stable_df = pd.DataFrame(
    sorted(feature_counts.items(), key=lambda x: (-x[1], x[0])),
    columns=["probe_id", "count"]
)

stable_df.to_csv("data/processed/stable_probe_counts.csv", index=False)


gpl = pd.read_csv(
    "data/GPL570-tbl-1.txt",
    sep="\t",
    comment="!",
    header=None,
    names=gpl_cols,
    low_memory=False,
)

gpl_small = gpl[["ID", "Gene Symbol", "Gene Title"]]
gpl_small = gpl_small.rename(columns={"ID": "probe_id"})

annotated = stable_df.merge(gpl_small, on="probe_id", how="left")

stable_mapped = annotated[
    (annotated["count"] >= 3) &
    (annotated["Gene Symbol"].notna()) &
    (annotated["Gene Symbol"].astype(str).str.strip() != "")
].copy()

clf.fit(X, y)

model = clf.named_steps["model"]
coef = model.coef_.flatten()

var = clf.named_steps["var"]
mask = var.get_support()

selected_features = np.array(feature_cols)[mask]

coef_df = pd.DataFrame({
    "probe_id": selected_features,
    "coef": coef
})
coef_df["abs_coef"] = coef_df["coef"].abs()

ranked = stable_mapped.merge(coef_df, on="probe_id", how="left")

ranked = ranked.sort_values(
    ["count", "abs_coef"],
    ascending=[False, False]
)


signature = ranked.head(15)["probe_id"].tolist()

#print("Signature genes:", signature)


X_sig = data[signature]

scores_sig = cross_val_score(
    clf,
    X_sig,
    y,
    cv=5,
    scoring="roc_auc"
)

'''print("Signature AUC:", scores_sig.mean(), "+-", scores_sig.std())



print("Full model AUC:", scores.mean())
print("Signature AUC:", scores_sig.mean())


print("\n--- Signature size sweep ---")'''

results = []

for k in range(2, 15):  # you can extend later
    sig = ranked.head(k)["probe_id"].tolist()
    X_sig = data[sig]

    scores_k = cross_val_score(
        clf,
        X_sig,
        y,
        cv=5,
        scoring="roc_auc"
    )

    mean_k = scores_k.mean()
    std_k = scores_k.std()

    results.append((k, mean_k, std_k))

    print(f"{k:2d} genes -> AUC: {mean_k:.3f} ± {std_k:.3f}")