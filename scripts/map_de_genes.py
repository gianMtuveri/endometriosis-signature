import pandas as pd

# ---------- LOAD PROBES ----------
up = pd.read_csv("data/processed/up_genes.txt", header=None, names=["probe_id"])
down = pd.read_csv("data/processed/down_genes.txt", header=None, names=["probe_id"])

# ---------- LOAD GPL ANNOTATION ----------
gpl_cols = [
    "ID","GB_ACC","SPOT_ID","Species Scientific Name","Annotation Date",
    "Sequence Type","Sequence Source","Target Description",
    "Representative Public ID","Gene Title","Gene Symbol",
    "ENTREZ_GENE_ID","RefSeq Transcript ID",
    "Gene Ontology Biological Process",
    "Gene Ontology Cellular Component",
    "Gene Ontology Molecular Function",
]

gpl = pd.read_csv(
    "data/raw/GPL570-tbl-1.txt",
    sep="\t",
    comment="!",
    header=None,
    names=gpl_cols,
    low_memory=False,
)

gpl_small = gpl[["ID", "Gene Symbol"]].rename(columns={"ID": "probe_id"})

# ---------- MAP ----------
up_mapped = up.merge(gpl_small, on="probe_id", how="left")
down_mapped = down.merge(gpl_small, on="probe_id", how="left")

# ---------- CLEAN ----------
def clean(df):
    df = df[df["Gene Symbol"].notna()]
    df = df[df["Gene Symbol"].str.strip() != ""]
    df["Gene Symbol"] = df["Gene Symbol"].str.split(" /// ").str[0]
    return df.drop_duplicates(subset="Gene Symbol")

up_clean = clean(up_mapped)
down_clean = clean(down_mapped)

# ---------- SAVE ----------
up_clean["Gene Symbol"].to_csv("data/processed/up_genes_mapped.txt", index=False, header=False)
down_clean["Gene Symbol"].to_csv("data/processed/down_genes_mapped.txt", index=False, header=False)

print("UP genes:", len(up_clean))
print("DOWN genes:", len(down_clean))

print("\nTop UP genes:")
print(up_clean.head(10))

print("\nTop DOWN genes:")
print(down_clean.head(10))