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
from typing import Any

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


def read_parquet(file_path: Path | str) -> pd.DataFrame | pd.Series:
    """
    Read a parquet file into a pandas DataFrame or Series.

    Uses pyarrow engine with automatic index restoration.

    Parameters
    ----------
    file_path : Path | str
        Path to the parquet file

    Returns
    -------
    pd.DataFrame | pd.Series
        The loaded data with index preserved

    Raises
    ------
    FileNotFoundError
        If the file does not exist
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Parquet file not found: {file_path}")

    return pd.read_parquet(file_path, engine="pyarrow")


# =============================================================================
# WRITING
# =============================================================================
def write_parquet(
    data: pd.DataFrame | pd.Series,
    file_path: Path | str,
    compression: str = "snappy",
    **kwargs: Any,
) -> None:
    """
    Write a pandas DataFrame or Series to parquet format.

    Parameters
    ----------
    data : pd.DataFrame | pd.Series
        Data to write
    file_path : Path | str
        Output path (will be created with .parquet extension if missing)
    compression : str, default='snappy'
        Compression algorithm: 'snappy' (fast), 'gzip', 'brotli', 'zstd', or None
    **kwargs
        Additional arguments passed to to_parquet()

    Notes
    -----
    - Index is always preserved (index=True by default unless overridden in kwargs)
    - Uses pyarrow engine for compatibility and performance
    - Creates parent directories if needed
    - Snappy compression provides good balance of speed and compression ratio
    """
    file_path = Path(file_path)
    file_path.parent.mkdir(parents=True, exist_ok=True)

    # Ensure index is preserved by default
    if 'index' not in kwargs:
        kwargs['index'] = True

    data.to_parquet(file_path, engine="pyarrow", compression=compression, **kwargs)
    logger.debug(f"Wrote parquet: {file_path} ({file_path.stat().st_size / 1024:.1f} KB)")
