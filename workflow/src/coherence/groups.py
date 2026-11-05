"""Group selection shared by the coherence compute and null-cache stages.

The cache builder has to enumerate every group size (n_scored_members) that the
per-source runs will score, and it must do so with exactly the same filters the
compute stage applies. Keeping one copy of `build_groups` here is what makes that
guarantee real; a second, drifting copy would silently leave some n out of the
cache and make the compute stage fall back or fail.
"""

# =============================================================================
# IMPORTS
# =============================================================================
from __future__ import annotations

# 3. Third-party Imports
from loguru import logger

# 2. Data Processing Imports
import pandas as pd


# =============================================================================
# CORE LOGIC
# =============================================================================
def build_groups(
    background: pd.DataFrame,
    long_table: pd.DataFrame,
    min_group_size: int,
    max_group_size: int,
    max_term_genes: int,
) -> dict[str, pd.DataFrame]:
    """Map surviving groups (keyed on group_id) -> their DR<threshold member rows."""
    merged = long_table.merge(
        background[["Systematic ID", "norm_DR", "norm_DL"]], on="Systematic ID", how="inner"
    )
    groups = {}
    for group_id, grp in merged.groupby("group_id"):
        grp = grp.drop_duplicates(subset="Systematic ID")
        n_total = int(grp["n_annotated_members"].iloc[0])
        if n_total > max_term_genes:
            continue
        if min_group_size <= len(grp) <= max_group_size:
            groups[group_id] = grp
    logger.info(
        f"{merged['group_id'].nunique():,} groups with >=1 background member -> "
        f"{len(groups):,} with {min_group_size} <= size <= {max_group_size} "
        f"and n_annotated_members <= {max_term_genes}"
    )
    return groups
