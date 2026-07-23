"""Contract tests for coherence source adapters (unified long-table schema)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.src.coherence.sources import (
    LONG_TABLE_COLUMNS,
    load_macrocomplex,
)

LONG_TABLE_COLUMNS_EXPECTED = [
    "source", "group_id", "group_name", "Systematic ID", "Name", "n_group_genes",
]


def _write_macrocomplex(tmp_path: Path) -> Path:
    """A tiny macromolecular_complex_annotation.tsv: 2 complexes, one member with no symbol."""
    df = pd.DataFrame(
        {
            "complex_term_id": ["GO:0001", "GO:0001", "GO:0002", "GO:0002"],
            "GO_term_name": ["alpha complex", "alpha complex", "beta complex", "beta complex"],
            "systematic_id": ["SPAC1", "SPAC2", "SPBC1", "SPBC2"],
            "symbol": ["gene1", None, "gene3", "gene4"],
        }
    )
    d = tmp_path / "ontologies_and_associations"
    d.mkdir(parents=True)
    path = tmp_path / "ontologies_and_associations" / "macromolecular_complex_annotation.tsv"
    df.to_csv(path, sep="\t", index=False)
    return tmp_path


def test_long_table_columns_constant():
    assert list(LONG_TABLE_COLUMNS) == LONG_TABLE_COLUMNS_EXPECTED


def test_macrocomplex_returns_contract_columns(tmp_path):
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_macrocomplex(pombase_dir)
    assert list(out.columns) == LONG_TABLE_COLUMNS_EXPECTED
    assert set(out["source"]) == {"go_macrocomplex"}


def test_macrocomplex_fills_missing_name_with_systematic_id(tmp_path):
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_macrocomplex(pombase_dir)
    row = out[out["Systematic ID"] == "SPAC2"].iloc[0]
    assert row["Name"] == "SPAC2"  # no symbol -> filled with systematic id


def test_macrocomplex_n_group_genes_is_per_term_total(tmp_path):
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_macrocomplex(pombase_dir)
    alpha = out[out["group_name"] == "alpha complex"]
    assert (alpha["n_group_genes"] == 2).all()
