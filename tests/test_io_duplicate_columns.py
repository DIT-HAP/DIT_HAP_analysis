"""Tests for duplicate column name handling in write_parquet/read_parquet."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
from workflow.src.io_table import read_parquet, write_parquet


def test_write_and_read_duplicate_columns_preserves_them(tmp_path):
    """Duplicate column names are preserved through write_parquet → read_parquet round-trip."""
    # Create a DataFrame with duplicate column names (intentional byte-faithful quirk)
    df = pd.DataFrame(
        [[1, 4, 7], [2, 5, 8], [3, 6, 9]],
        columns=["A", "B", "A"],  # Duplicate column name
    )

    parquet_file = tmp_path / "duplicate_cols.parquet"
    write_parquet(df, parquet_file)

    df_read = read_parquet(parquet_file)

    # Verify duplicate column names are preserved
    assert df.columns.tolist() == ["A", "B", "A"]
    assert df_read.columns.tolist() == ["A", "B", "A"]
    assert not df.columns.is_unique
    assert not df_read.columns.is_unique

    # Verify values are preserved
    assert df.equals(df_read)


def test_duplicate_columns_preserve_distinct_values(tmp_path):
    """The two columns with the same name keep their distinct values."""
    df = pd.DataFrame(
        [[10, 30, 50], [20, 40, 60]],
        columns=["X", "Y", "X"],  # Second 'X' column has different values
    )

    parquet_file = tmp_path / "test.parquet"
    write_parquet(df, parquet_file)
    df_read = read_parquet(parquet_file)

    # Access by position: first 'X' is [10, 20], second 'X' is [50, 60].
    assert df.iloc[:, 0].tolist() == [10, 20]
    assert df.iloc[:, 1].tolist() == [30, 40]
    assert df.iloc[:, 2].tolist() == [50, 60]
    assert df_read.iloc[:, 0].tolist() == [10, 20]
    assert df_read.iloc[:, 1].tolist() == [30, 40]
    assert df_read.iloc[:, 2].tolist() == [50, 60]


def test_phenotype_level_duplicate_essentiality_column(tmp_path):
    """Phenotype-level intentional duplicate DeletionLibrary_essentiality is preserved."""
    # Simulate the phenotype-level structure with duplicate DeletionLibrary_essentiality
    df = pd.DataFrame(
        [
            ["viable", "V", "WT-like", "WT-like", 11, "V", "V"],
            ["viable", "E", "microcolonies", "microcolonies", 5, "E", "E"],
        ],
        columns=[
            "FYPOviability",
            "DeletionLibrary_essentiality",
            "DeletionLibrary_category",
            "Sub_category",
            "Growth_tier",
            "DeletionLibrary_essentiality",  # Intentional duplicate
            "RevisedDeletionLibrary_essentiality",
        ],
    )

    parquet_file = tmp_path / "phenotype.parquet"
    write_parquet(df, parquet_file)
    df_read = read_parquet(parquet_file)

    # Verify the duplicate column is preserved
    assert df.columns.tolist().count("DeletionLibrary_essentiality") == 2
    assert df_read.columns.tolist().count("DeletionLibrary_essentiality") == 2

    # Verify both copies have the same values (they're derived from the same source)
    assert df.iloc[:, 1].tolist() == df.iloc[:, 5].tolist()
    assert df_read.iloc[:, 1].tolist() == df_read.iloc[:, 5].tolist()


def test_no_duplicate_columns_works_normally(tmp_path):
    """DataFrames without duplicate columns continue to work as before."""
    df = pd.DataFrame(
        {
            "A": [1, 2],
            "B": [3, 4],
            "C": [5, 6],
        }
    )

    parquet_file = tmp_path / "normal.parquet"
    write_parquet(df, parquet_file)
    df_read = read_parquet(parquet_file)

    assert df.columns.is_unique
    assert df_read.columns.is_unique
    assert df.equals(df_read)
