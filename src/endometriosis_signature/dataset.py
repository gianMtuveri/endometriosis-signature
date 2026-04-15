from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

from .constants import CLEAN_SOURCES, META_COLUMNS


def load_expression_matrix(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="!")
    df = df.set_index(df.columns[0]).T
    df.index.name = "sample_id"
    df.columns.name = "gene"
    return df


def parse_family_xml(path: str | Path) -> pd.DataFrame:
    ns = {"geo": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}
    tree = ET.parse(path)
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
    meta["target"] = (meta["endo_status"] == "Endometriosis").astype(int)
    return meta


def merge_expression_and_metadata(
    expression_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
) -> pd.DataFrame:
    meta = metadata_df.set_index("sample_id")
    data = expression_df.join(
        meta[["target", "source", "tissue", "severity"]],
        how="inner",
    )
    return data


def get_feature_columns(data: pd.DataFrame) -> list[str]:
    return [c for c in data.columns if c not in META_COLUMNS]


def make_clean_subset(data: pd.DataFrame) -> pd.DataFrame:
    return data[data["source"].isin(CLEAN_SOURCES)].copy()


def save_feature_columns(feature_cols: list[str], path: str | Path) -> None:
    with open(path, "w") as f:
        json.dump(feature_cols, f)


def save_dataset(data: pd.DataFrame, path: str | Path) -> None:
    data.to_parquet(path)