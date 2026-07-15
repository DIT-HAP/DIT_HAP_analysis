"""
Generic Table Readers
======================

File-extension-dispatched table loading, factored out of the original
`workflow/src/utils.py` (DIT_HAP_pipeline). No biology-specific logic —
safe to import from any module in this repository.

Input
-----
- Any .tsv/.bed/.csv/.xlsx file path

Output
------
- pandas DataFrame

Usage
-----
    from workflow.src.io import read_file
    df = read_file(Path("resources/curated/essentiality_verification.csv"))

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# =============================================================================
# CORE LOGIC
# =============================================================================
def read_file(file: Path, **kwargs) -> pd.DataFrame:
    """Read a table into a DataFrame, dispatching on file extension (tsv/bed/csv/xlsx)."""
    if "tsv" in file.name:
        return pd.read_csv(file, sep="\t", **kwargs)
    elif "bed" in file.name:
        return pd.read_csv(file, sep="\t", **kwargs)
    elif "csv" in file.name:
        return pd.read_csv(file, sep=",", **kwargs)
    elif "xlsx" in file.name:
        return pd.read_excel(file, **kwargs)
    else:
        raise ValueError(f"Unsupported file type: {file.name}")
