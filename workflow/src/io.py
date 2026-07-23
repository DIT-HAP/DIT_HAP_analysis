"""
Generic Table Readers and Writers
==================================

File-extension-dispatched table loading and saving, factored out of the original
`workflow/src/utils.py` (DIT_HAP_pipeline). Includes parquet support for efficient
intermediate file storage. No biology-specific logic — safe to import from any
module in this repository.

Input
-----
- Any .tsv/.bed/.csv/.xlsx/.parquet file path

Output
------
- pandas DataFrame or Series

Usage
-----
    from workflow.src.io import read_file, read_parquet, write_parquet
    df = read_file(Path("resources/curated/essentiality_verification.csv"))
    write_parquet(df, Path("results/intermediate/data.parquet"))
    df2 = read_parquet(Path("results/intermediate/data.parquet"))

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# =============================================================================
# READING
# =============================================================================
def read_file(file: Path, **kwargs) -> pd.DataFrame:
    """Read a table into a DataFrame, dispatching on file extension (tsv/bed/csv/xlsx/parquet)."""
    if "tsv" in file.name:
        return pd.read_csv(file, sep="\t", **kwargs)
    elif "bed" in file.name:
        return pd.read_csv(file, sep="\t", **kwargs)
    elif "csv" in file.name:
        return pd.read_csv(file, sep=",", **kwargs)
    elif "xlsx" in file.name:
        return pd.read_excel(file, **kwargs)
    elif "parquet" in file.name:
        return pd.read_parquet(file, engine="pyarrow", **kwargs)
    else:
        raise ValueError(f"Unsupported file type: {file.name}")


# Schema-metadata keys used to make the Series <-> DataFrame round-trip transparent.
# Parquet only stores tables, and pd.Series has no .to_parquet(), so write_parquet
# frames a Series and stamps these markers; read_parquet reads them to squeeze back.
_OBJ_TYPE_KEY = b"__pandas_object_type__"
_SERIES_HAS_NAME_KEY = b"__series_has_name__"
_SERIES_PLACEHOLDER_COL = "__series_values__"


def read_parquet(file_path: Path | str) -> pd.DataFrame | pd.Series:
    """
    Read a parquet file into a pandas DataFrame or Series.

    Uses pyarrow with automatic index restoration. If the file was written from a
    pd.Series by write_parquet (marked in the schema metadata), it is squeezed back
    to a Series with its original name; otherwise a DataFrame is returned.

    Parameters
    ----------
    file_path : Path | str
        Path to the parquet file

    Returns
    -------
    pd.DataFrame | pd.Series
        The loaded data with index (and, for a Series, name) preserved

    Raises
    ------
    FileNotFoundError
        If the file does not exist
    """
    import pyarrow.parquet as pq

    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Parquet file not found: {file_path}")

    table = pq.read_table(file_path)
    metadata = table.schema.metadata or {}
    obj_type = metadata.get(_OBJ_TYPE_KEY, b"dataframe")

    frame = table.to_pandas()

    # Restore duplicate column names if they were renamed during write
    if b"__original_columns__" in metadata:
        import ast
        original_columns = ast.literal_eval(metadata[b"__original_columns__"].decode("utf-8"))
        frame.columns = original_columns

    if obj_type != b"series":
        return frame

    # Squeeze the single stored column back into a Series, restoring its name.
    series = frame.iloc[:, 0]
    if metadata.get(_SERIES_HAS_NAME_KEY, b"0") == b"1":
        series.name = frame.columns[0]
    else:
        series.name = None
    return series


# =============================================================================
# WRITING
# =============================================================================
def write_parquet(
    data: pd.DataFrame | pd.Series,
    file_path: Path | str,
    compression: str = "snappy",
) -> None:
    """
    Write a pandas DataFrame or Series to parquet format.

    Parameters
    ----------
    data : pd.DataFrame | pd.Series
        Data to write. A Series is stored as a one-column table and stamped in the
        schema metadata so read_parquet can squeeze it back to a Series (with its
        original name). pd.Series has no native .to_parquet(), hence this framing.
    file_path : Path | str
        Output path
    compression : str, default='snappy'
        Compression algorithm: 'snappy' (fast), 'gzip', 'brotli', 'zstd', or None

    Notes
    -----
    - The index is always preserved (stored via pyarrow's preserve_index=True)
    - Uses pyarrow for compatibility and performance
    - Creates parent directories if needed
    - Snappy compression provides a good balance of speed and compression ratio
    - Handles duplicate column names by suffixing them for storage, then restoring
      on read (intentional byte-faithful quirk for phenotype features)
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    file_path = Path(file_path)
    file_path.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(data, pd.Series):
        has_name = data.name is not None
        col_name = data.name if has_name else _SERIES_PLACEHOLDER_COL
        frame = data.to_frame(name=col_name)
        extra_metadata = {
            _OBJ_TYPE_KEY: b"series",
            _SERIES_HAS_NAME_KEY: b"1" if has_name else b"0",
        }
    else:
        frame = data
        extra_metadata = {_OBJ_TYPE_KEY: b"dataframe"}

    # Handle duplicate column names: pyarrow rejects them, but our phenotype-level
    # features intentionally carry a duplicate "DeletionLibrary_essentiality" column
    # to match the reference output byte-for-byte. Rename duplicates for storage,
    # stamp the original names in metadata, and restore them on read.
    if isinstance(frame, pd.DataFrame) and not frame.columns.is_unique:
        original_columns = frame.columns.tolist()
        # Suffix duplicates: DeletionLibrary_essentiality → DeletionLibrary_essentiality__dup1, etc.
        new_columns = []
        seen = {}
        for col in original_columns:
            if col in seen:
                seen[col] += 1
                new_columns.append(f"{col}__dup{seen[col]}")
            else:
                seen[col] = 0
                new_columns.append(col)
        frame = frame.copy()
        frame.columns = new_columns
        extra_metadata[b"__original_columns__"] = str(original_columns).encode("utf-8")
    else:
        extra_metadata = {_OBJ_TYPE_KEY: b"dataframe"}

    table = pa.Table.from_pandas(frame, preserve_index=True)
    # Merge our markers into the existing (pandas) schema metadata rather than replace.
    merged_metadata = {**(table.schema.metadata or {}), **extra_metadata}
    table = table.replace_schema_metadata(merged_metadata)

    pq.write_table(table, file_path, compression=compression)
    logger.debug(f"Wrote parquet: {file_path} ({file_path.stat().st_size / 1024:.1f} KB)")
