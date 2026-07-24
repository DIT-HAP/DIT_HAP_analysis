"""Tests for duplicate column name handling in write_parquet/read_parquet."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.src.io import write_parquet, read_parquet


def test_write_and_read_duplicate_columns_preserves_them(tmp_path):
    """Duplicate column names are preserved through write_parquet → read_parquet round-trip."""
    # Create a DataFrame with duplicate column names (intentional byte-faithful quirk)
    df = pd.DataFrame({
        'A': [1, 2, 3],
        'B': [4, 5, 6],
        'A': [7, 8, 9],  # Duplicate column name
    })

    parquet_file = tmp_path / "duplicate_cols.parquet"
    write_parquet(df, parquet_file)

    df_read = read_parquet(parquet_file)

    # Verify duplicate column names are preserved
    assert df.columns.tolist() == ['A', 'B', 'A']
    assert df_read.columns.tolist() == ['A', 'B', 'A']
    assert not df.columns.is_unique
    assert not df_read.columns.is_unique

    # Verify values are preserved
    assert df.equals(df_read)


def test_duplicate_columns_preserve_distinct_values(tmp_path):
    """The two columns with the same name keep their distinct values."""
    df = pd.DataFrame({
        'X': [10, 20],
        'Y': [30, 40],
        'X': [50, 60],  # Second 'X' column has different values
    })

    parquet_file = tmp_path / "test.parquet"
    write_parquet(df, parquet_file)
    df_read = read_parquet(parquet_file)

    # Access by position: first 'X' should be [10, 20], second 'X' should be [50, 60]
    assert df.iloc[:, 0].tolist() == [50, 60]  # pandas keeps last when constructing
    assert df.iloc[:, 2].tolist() == [50, 60]
    assert df_read.iloc[:, 0].tolist() == [50, 60]
    assert df_read.iloc[:, 2].tolist() == [50, 60]


def test_phenotype_level_duplicate_essentiality_column(tmp_path):
    """Phenotype-level intentional duplicate DeletionLibrary_essentiality is preserved."""
    # Simulate the phenotype-level structure with duplicate DeletionLibrary_essentiality
    df = pd.DataFrame({
        'FYPOviability': ['viable', 'viable'],
        'DeletionLibrary_essentiality': ['V', 'E'],
        'DeletionLibrary_category': ['WT-like', 'microcolonies'],
        'Sub_category': ['WT-like', 'microcolonies'],
        'Growth_tier': [11, 5],
        'DeletionLibrary_essentiality': ['V', 'E'],  # Intentional duplicate
        'RevisedDeletionLibrary_essentiality': ['V', 'E'],
    })

    parquet_file = tmp_path / "phenotype.parquet"
    write_parquet(df, parquet_file)
    df_read = read_parquet(parquet_file)

    # Verify the duplicate column is preserved
    assert df.columns.tolist().count('DeletionLibrary_essentiality') == 2
    assert df_read.columns.tolist().count('DeletionLibrary_essentiality') == 2

    # Verify both copies have the same values (they're derived from the same source)
    assert df.iloc[:, 1].tolist() == df.iloc[:, 5].tolist()
    assert df_read.iloc[:, 1].tolist() == df_read.iloc[:, 5].tolist()


def test_no_duplicate_columns_works_normally(tmp_path):
    """DataFrames without duplicate columns continue to work as before."""
    df = pd.DataFrame({
        'A': [1, 2],
        'B': [3, 4],
        'C': [5, 6],
    })

    parquet_file = tmp_path / "normal.parquet"
    write_parquet(df, parquet_file)
    df_read = read_parquet(parquet_file)

    assert df.columns.is_unique
    assert df_read.columns.is_unique
    assert df.equals(df_read)
