import pandas as pd
from scipy.stats import ttest_ind

data = pd.read_parquet("data/processed/endometriosis_clean_filtered.parquet")

meta_cols = ["target", "source", "tissue", "severity"]
feature_cols = [c for c in data.columns if c not in meta_cols]

X = data[feature_cols]
y = data["target"]

# split
X_case = X[y == 1]
X_control = X[y == 0]

rows = []

for gene in feature_cols:
    case_vals = X_case[gene]
    ctrl_vals = X_control[gene]

    stat, pval = ttest_ind(case_vals, ctrl_vals, equal_var=False)

    logfc = case_vals.mean() - ctrl_vals.mean()

    rows.append({
        "gene": gene,
        "logFC": logfc,
        "pval": pval
    })

df = pd.DataFrame(rows)
df["abs_logFC"] = df["logFC"].abs()

df = df.sort_values("pval")

df.to_csv("data/processed/differential_expression.csv", index=False)

print(df.head(20))


up = df[df["logFC"] > 0].head(100)
down = df[df["logFC"] < 0].head(100)


up["gene"].to_csv("data/processed/up_genes.txt", index=False, header=False)
down["gene"].to_csv("data/processed/down_genes.txt", index=False, header=False)

