"""
PCR / Library-Prep Quality Control — Core Logic
================================================

Shared loader for the PCR QC stage. Ported from
DIT_HAP_pipeline/workflow/notebooks/thesis_figures.ipynb ("2. PCR quality
control") and factored out of the original single-script port so the stage
can be split into independent Snakemake rules (prepare -> plot), each
re-runnable on its own.

Figure panels live in figure_render/: panels (a)-(c) use the generic scatter
renderer, panel (d) the spike-in renderer (figure_render/spikein.py).

Usage
-----
    from pcr_qc.core import read_merged_reads
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd


# =============================================================================
# LOADERS
# =============================================================================
def read_merged_reads(path: Path) -> pd.DataFrame:
    """Read a merged reads TSV indexed by (Chr, Coordinate, Strand) with PBL/PBR/Reads."""
    return pd.read_csv(path, sep="\t", index_col=[0, 1, 2])
