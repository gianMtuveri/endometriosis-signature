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

data = pd.read_parquet("data/processed/endometriosis_clean_filtered.parquet")

meta_cols = ["target", "source", "tissue", "severity"]
feature_cols = [col for col in data.columns if col not in meta_cols]

X = data[feature_cols]
y = data["target"]


print(X.shape, y.shape)

selector = VarianceThreshold(threshold=0.01)
X_sel = selector.fit_transform(X)



clf = Pipeline([
    ("var", VarianceThreshold(threshold=0.01)),
    ("scaler", StandardScaler()),
    ("model", LogisticRegression(max_iter=2000))
])

scores = cross_val_score(
    clf,
    X,
    y,
    cv=5,
    scoring="roc_auc"
)

print("CV AUC:", scores)
print("Mean:", scores.mean())
print("Std:", scores.std())

########## feature importance #########

clf.fit(X, y)

model = clf.named_steps["model"]
coef = model.coef_.flatten()

var = clf.named_steps["var"]
mask = var.get_support()

selected_features = np.array(feature_cols)[mask]

top_idx = np.argsort(np.abs(coef))[-20:]
top_idx = top_idx[np.argsort(np.abs(coef[top_idx]))[::-1]]

top_features = selected_features[top_idx]
top_weights = coef[top_idx]

'''print("\nTop features from full-data fit:")
for f, w in zip(top_features, top_weights):
    print(f, w)'''

########## feature stability across CV folds #########

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

'''print("\nTop stable features across folds:")
for feat, count in feature_counts.most_common(20):
    print(feat, count)'''

stable_feats = [feat for feat, count in feature_counts.items() if count >= 4]
#print(stable_feats)

stable_df = pd.DataFrame(
    sorted(feature_counts.items(), key=lambda x: (-x[1], x[0])),
    columns=["probe_id", "count"]
)
stable_df.to_csv("data/processed/stable_probe_counts.csv", index=False)


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

gpl = pd.read_csv(
    "data/GPL570-tbl-1.txt",
    sep="\t",
    comment="!",
    header=None,
    names=gpl_cols,
    low_memory=False,
)

#print(gpl.head())
#print(gpl[["ID", "Gene Symbol", "Gene Title"]].head())

gpl_small = gpl[["ID", "Gene Symbol", "Gene Title"]].copy()
gpl_small = gpl_small.rename(columns={"ID": "probe_id"})

stable_df = pd.read_csv("data/processed/stable_probe_counts.csv")

annotated = stable_df.merge(gpl_small, on="probe_id", how="left")
top = annotated.sort_values(["count", "probe_id"], ascending=[False, True]).head(20)

#print(top[["probe_id", "count", "Gene Symbol", "Gene Title"]])

stable_mapped = annotated[
    (annotated["count"] >= 3) &
    (annotated["Gene Symbol"].notna()) &
    (annotated["Gene Symbol"].astype(str).str.strip() != "")
].copy()

stable_mapped = stable_mapped.sort_values(
    ["count", "probe_id"],
    ascending=[False, True]
)

#print(stable_mapped[["probe_id", "count", "Gene Symbol", "Gene Title"]])


stable_mapped.to_csv("data/processed/stable_mapped_probes.csv", index=False)


coef_df = pd.DataFrame({
    "probe_id": selected_features,
    "coef": coef
})
coef_df["abs_coef"] = coef_df["coef"].abs()

ranked = stable_df.merge(coef_df, on="probe_id", how="left")
ranked = ranked.merge(gpl_small, on="probe_id", how="left")

ranked = ranked.sort_values(
    ["count", "abs_coef"],
    ascending=[False, False]
)

print(ranked[["probe_id", "count", "coef", "Gene Symbol", "Gene Title"]].head(20))