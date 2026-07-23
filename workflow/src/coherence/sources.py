"""Source adapters mapping PomBase grouping databases to a unified coherence long-table.

Every adapter returns the same contract columns (LONG_TABLE_COLUMNS), so the
downstream compute/plot stages are source-agnostic. Add a database = add one
adapter here + one entry in SOURCE_LOADERS + one line in config.coherence.sources.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

LONG_TABLE_COLUMNS = [
    "source", "group_id", "group_name", "Systematic ID", "Name", "n_group_genes",
]

# PomBase macrocomplex annotation column -> canonical contract name.
_MACRO_RENAME = {
    "complex_term_id": "group_id",
    "GO_term_name": "group_name",
    "systematic_id": "Systematic ID",
    "symbol": "Name",
}


def _finalize(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Fill missing Name with Systematic ID, add source + n_group_genes, order columns."""
    df = df.copy()
    df["source"] = source
    df["Name"] = df["Name"].fillna(df["Systematic ID"])
    # Dedup + count on group_id (the stable GO term ID), NOT group_name. The old
    # compute_complex_coherence.py grouped on GO_term_name; keying on the ID is safer
    # (two distinct term IDs could share a name) and matches how the compute stage
    # later forms groups, keeping n_group_genes consistent with the DR-member count.
    df = df.drop_duplicates(subset=["group_id", "Systematic ID"])
    counts = df.groupby("group_id")["Systematic ID"].transform("size")
    df["n_group_genes"] = counts
    return df[LONG_TABLE_COLUMNS].reset_index(drop=True)


def load_macrocomplex(pombase_dir: Path) -> pd.DataFrame:
    """Flat PomBase macromolecular_complex_annotation.tsv -> unified long-table."""
    path = Path(pombase_dir) / "ontologies_and_associations" / "macromolecular_complex_annotation.tsv"
    raw = pd.read_csv(path, sep="\t").rename(columns=_MACRO_RENAME)
    for required in ["group_id", "group_name", "Systematic ID"]:
        if required not in raw.columns:
            raise ValueError(f"macrocomplex annotation missing '{required}' (have: {list(raw.columns)})")
    if "Name" not in raw.columns:
        raw["Name"] = pd.NA
    return _finalize(raw[["group_id", "group_name", "Systematic ID", "Name"]], "go_macrocomplex")
