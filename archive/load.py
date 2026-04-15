import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import cross_val_score


df = pd.read_csv("data/GSE51981_series_matrix.txt", sep="\t", comment="!")
#print(df.shape)
#print(df.head())
#print(df.columns.tolist())

df = df.set_index(df.columns[0])  # set gene column as index
df = df.T  # transpose

'''print(df.head())
print(df.shape)
print(df.isna().sum().sum())
print(df.dtypes.unique())     # should be floatsum())

print(df.iloc[:5, :5])
df.index.name = "sample_id"
df.columns.name = "gene"'''



xml_path = "data/GSE51981_family.xml"
ns = {"geo": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}

tree = ET.parse(xml_path)
root = tree.getroot()

rows = []

for sample in root.findall("geo:Sample", ns):
    sample_id = sample.get("iid")
    title = sample.findtext("geo:Title", default=None, namespaces=ns)
    accession = sample.findtext("geo:Accession", default=None, namespaces=ns)
    description = sample.findtext("geo:Description", default=None, namespaces=ns)

    channel = sample.find("geo:Channel", ns)

    source = None
    tissue = None
    endo_status = None
    severity = None

    if channel is not None:
        source = channel.findtext("geo:Source", default=None, namespaces=ns)

        for char in channel.findall("geo:Characteristics", ns):
            tag = char.get("tag")
            value = (char.text or "").strip()

            if tag == "tissue":
                tissue = value
            elif tag == "endometriosis/no endometriosis":
                endo_status = value
            elif tag == "endometriosis severity":
                severity = value

    rows.append(
        {
            "sample_id": sample_id,
            "title": title,
            "accession": accession,
            "description": description,
            "source": source,
            "tissue": tissue,
            "endo_status": endo_status,
            "severity": severity,
        }
    )

meta = pd.DataFrame(rows)

#print(meta.head())
#print(meta["source"].value_counts(dropna=False))
#print(meta["endo_status"].value_counts(dropna=False))
#print(meta["severity"].value_counts(dropna=False))

meta["target"] = (meta["endo_status"] == "Endometriosis").astype(int)
#print(meta["target"].value_counts())



print(df.index[:5])
print(meta["sample_id"].head())

meta = meta.set_index("sample_id")


common = df.index.intersection(meta.index)
print(len(common))
print(len(df), len(meta))


data = df.join(meta[["target", "source", "tissue", "severity"]], how="inner")

data.to_parquet("data/processed/endometriosis_clean.parquet")

print(data.shape)
print(data[["target"]].head())
print(data["target"].value_counts())

missing = set(df.index) - set(data.index)
print(len(missing))


print(data.groupby("target")["source"].value_counts())



feature_cols = df.columns  # genes only

X = data[feature_cols]
y = data["target"]

X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size=0.2,
    stratify=y,
    random_state=42
)

clf = Pipeline([
    ("scaler", StandardScaler()),
    ("model", LogisticRegression(max_iter=2000))
])

clf.fit(X_train, y_train)

proba = clf.predict_proba(X_test)[:, 1]

print("ROC-AUC:", roc_auc_score(y_test, proba))

print(data.groupby("target")["tissue"].value_counts(normalize=True))
print(data.groupby("target")["source"].value_counts(normalize=True))



clean_data = data[
    data["source"].isin([
        "Endometriosis_Moderate/Severe",
        "Endometriosis_Minimal/Mild",
        "Non-Endometriosis_No Uterine Pelvic Pathology"
    ])
]

clean_data.to_parquet("data/processed/endometriosis_clean_filtered.parquet")

print(clean_data["target"].value_counts())

X_clean = clean_data[feature_cols]
y_clean = clean_data["target"]

from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X_clean, y_clean,
    test_size=0.2,
    stratify=y_clean,
    random_state=42
)

clf.fit(X_train, y_train)

proba = clf.predict_proba(X_test)[:, 1]

from sklearn.metrics import roc_auc_score
print("ROC-AUC clean:", roc_auc_score(y_test, proba))


train_proba = clf.predict_proba(X_train)[:, 1]
test_proba = clf.predict_proba(X_test)[:, 1]

print("Train AUC:", roc_auc_score(y_train, train_proba))
print("Test AUC:", roc_auc_score(y_test, test_proba))


scores = cross_val_score(
    clf,
    X_clean,
    y_clean,
    cv=5,
    scoring="roc_auc"
)

print("CV AUC:", scores)
print("Mean:", scores.mean())

y_shuffled = np.random.permutation(y_clean)

scores_shuff = cross_val_score(
    clf,
    X_clean,
    y_shuffled,
    cv=5,
    scoring="roc_auc"
)

print("Shuffled AUC:", scores_shuff.mean())

X_clean = clean_data[feature_cols]
y_clean = clean_data["target"]

from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X_clean, y_clean,
    test_size=0.2,
    stratify=y_clean,
    random_state=42
)

clf.fit(X_train, y_train)

proba = clf.predict_proba(X_test)[:, 1]

from sklearn.metrics import roc_auc_score
print("ROC-AUC clean:", roc_auc_score(y_test, proba))

with open("data/processed/feature_columns.json", "w") as f:
    json.dump(feature_cols, f)