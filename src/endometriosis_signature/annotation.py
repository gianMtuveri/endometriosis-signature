from __future__ import annotations

from pathlib import Path

import pandas as pd

from .constants import GPL570_COLUMNS


def load_gpl570_annotation(path: str | Path) -> pd.DataFrame:
    gpl = pd.read_csv(
        path,
        sep="\t",
        comment="!",
        header=None,
        names=GPL570_COLUMNS,
        low_memory=False,
    )
    gpl_small = gpl[["ID", "Gene Symbol", "Gene Title"]].copy()
    gpl_small = gpl_small.rename(columns={"ID": "probe_id"})
    return gpl_small


def annotate_probes(stable_df: pd.DataFrame, gpl_small: pd.DataFrame) -> pd.DataFrame:
    return stable_df.merge(gpl_small, on="probe_id", how="left")