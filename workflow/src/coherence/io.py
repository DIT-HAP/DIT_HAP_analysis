"""Shared input loaders for the coherence chain.

Every coherence stage reads the same two upstream tables, and each used to carry
its own copy of the loading rules — the um/lam rename, the `Systematic ID`
first-column normalization, the +/-inf scrub, the DR threshold, and the
`(DR, DL/10)` divisors. The divisors additionally had a second copy in
`coherence.metrics`. They live here now, so one schema or normalization change
has one place to land.

`LONG_TABLE_COLUMNS` is re-exported from `coherence.sources` rather than
re-declared: that module is what *produces* the long-table, so it owns the
contract.

Input
-----
- Nothing: a contract re-export plus two frame-in/frame-out loaders.

Output
------
- `FITTING_RESULTS_COLUMNS`, `LONG_TABLE_COLUMNS`
- `load_fitting_results()`, `load_long_table()`
"""
# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from __future__ import annotations

from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
from coherence.metrics import DL_NORM_MAX, DR_NORM_MAX
from coherence.sources import LONG_TABLE_COLUMNS
from release_schema import normalize_legacy_metrics


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# The upstream per-gene fitness table's required columns, after the index column
# has been normalized to `Systematic ID` and um/lam have been renamed to DR/DL.
FITTING_RESULTS_COLUMNS = ["Systematic ID", "DR", "DL"]


# =============================================================================
# CORE LOGIC
# =============================================================================
def load_fitting_results(path: Path, *, dr_threshold: float | None = None) -> pd.DataFrame:
    """Load upstream fitting_results.tsv into per-gene DR/DL fitness points."""
    # Handles the three upstream quirks in one place: the systematic id arrives as
    # the first column under whatever name the release used (the index here; see
    # `clustering/candidates.py`, which reads the same table), legacy releases still
    # ship the pre-rename `um`/`lam` headers, and DR/DL may be non-finite.
    # `path` is the upstream `gene_level/fitting_results.tsv`.
    #
    # norm_DR / norm_DL are DR / DR_NORM_MAX and DL / DL_NORM_MAX. Those divisors
    # are plain constants, not a fitted range, so this is an affine reshape of
    # (DR, DL/10) — the coherence metrics are all Euclidean distances in that
    # space, which a reflection of one axis leaves unchanged. Upstream flipped the
    # DR sign on 2026-09-17 (negative = depleted), mirroring the space without
    # moving any z-score; only which side of the DR axis WT sits on changes.
    #
    # dr_threshold keeps only genes with DR < dr_threshold (negative = depleted,
    # matching the upstream sign convention) and logs the reduction.
    fitting = pd.read_csv(path, sep="\t", index_col=0).reset_index()
    if "Systematic ID" not in fitting.columns:
        first_col = fitting.columns[0]
        logger.info(f"Renaming fitting_results index column '{first_col}' -> 'Systematic ID'")
        fitting = fitting.rename(columns={first_col: "Systematic ID"})

    fitting = normalize_legacy_metrics(fitting)
    missing = [col for col in FITTING_RESULTS_COLUMNS if col not in fitting.columns]
    if missing:
        raise ValueError(
            f"fitting_results.tsv missing required column(s) {missing} (have: {list(fitting.columns)})"
        )

    # Map +/-inf to NaN before dropna so non-finite DR/DL never reach a distance
    # computation (dropna alone keeps +/-inf).
    fitting = fitting.replace([np.inf, -np.inf], np.nan).dropna(subset=["DR", "DL"]).copy()
    fitting["norm_DR"] = fitting["DR"].to_numpy(dtype=float) / DR_NORM_MAX
    fitting["norm_DL"] = fitting["DL"].to_numpy(dtype=float) / DL_NORM_MAX

    if dr_threshold is None:
        return fitting

    background = fitting[fitting["DR"] < dr_threshold].copy()
    logger.info(
        f"fitting_results.tsv: {len(fitting):,} fitted genes -> "
        f"{len(background):,} background genes with DR < {dr_threshold}"
    )
    return background


def load_long_table(path: Path) -> pd.DataFrame:
    """Load a prepared group_annotation_long.tsv and check its contract."""
    # `workflow/src/coherence/sources.py` is the only producer, so the contract
    # lives there; this only asserts it, once, for every consumer.
    long_table = pd.read_csv(path, sep="\t")
    missing = [col for col in LONG_TABLE_COLUMNS if col not in long_table.columns]
    if missing:
        raise ValueError(
            f"annotation long-table missing required column(s) {missing} (have: {list(long_table.columns)})"
        )
    return long_table
