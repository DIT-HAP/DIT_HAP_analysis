"""Per-group fractions the coherence tables and figures share (theme D).

Two descriptive quantities a term carries, computed here so a table and the
figure drawn from it can never disagree:

(a) `paralog_fraction` — the share of a group's members whose feature-matrix
    `paralog_count` is > 0. A high value flags a group whose members' deletion
    phenotypes may be buffered by a redundant paralog (DR dampened toward WT),
    which is a candidate explanation for a group looking internally dispersed.
(b) `moonlighting_fraction` — the share of a term's members that sit in MORE
    de-duplication groups than the mode. After de-duplication a high group count
    is not GO nesting: the gene genuinely keeps turning up in modules that are
    not redundant with each other, so a term carried by such genes may be
    incoherent because several pathways pull on it at once.

Breadth is per VIEW (one source, `combined`, or `dedup`). A view's groups are the
de-duplication clusters its terms belong to, and a gene's count is how many of
those clusters contain it. The mode of that per-view count distribution is the
view's own cut, so every view is read against its own breadth profile rather than
against one global number — the gene sets themselves stay the cluster's own, so
each view counts the same groups the same way and only the selection differs.

Pure functions over frames; no IO.
"""
# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from __future__ import annotations

from collections.abc import Iterable

# 2. Data Processing Imports
import numpy as np
import pandas as pd


# =============================================================================
# CORE LOGIC — paralog carry-over
# =============================================================================
def paralog_ids_from_features(features: pd.DataFrame) -> set[str]:
    """Genes with >=1 paralog, from the feature matrix's `paralog_count` column."""
    # Whichever source fills that column is `features.paralog_source` in
    # config/analysis.yaml (the deletion library, currently); nothing here picks a
    # source of its own, so this set and the ML feature cannot drift apart. A gene
    # absent from the matrix (not a coding gene, or outside the deletion library)
    # reads as "no paralog", which is what the column itself says.
    return set(features.loc[features["paralog_count"] > 0, "gene_systematic_id"].astype(str))


def paralog_fraction(members: Iterable[str], paralog_ids: set[str]) -> float:
    """Fraction of `members` that have >=1 paralog (present in paralog_ids)."""
    members = list(members)
    if not members:
        return np.nan
    return sum(1 for member in members if member in paralog_ids) / len(members)


# =============================================================================
# CORE LOGIC — gene breadth (moonlighting)
# =============================================================================
def view_breadth(group_members: pd.DataFrame, clusters: set[str]) -> tuple[pd.DataFrame, int]:
    """Per-gene group count within `clusters`, plus the mode that cuts it.

    `group_members` is the de-duplication stage's (group, gene) long table: one row
    per group per gene, carrying `redundancy_cluster`. Restricting it to a view's
    clusters is the whole of "this view's breadth" — the gene sets stay each
    cluster's own, so two views differ only in which groups they count.

    Returns `(per_gene, mode)` with per_gene = gene / n_groups / is_moonlighting,
    ordered by breadth descending so the broadest genes read first.
    """
    view = group_members[group_members["redundancy_cluster"].isin(clusters)]
    # One row per group per gene, so a gene's row count IS its group count.
    counts = view.groupby("gene").size().rename("n_groups").reset_index()
    # The mode is over GENES, not rows: a gene's row count grows with its own group
    # count, so a row-weighted mode lands on the wrong value and the cut it defines
    # stops meaning "more groups than the typical gene".
    mode = int(counts["n_groups"].mode().iloc[0]) if not counts.empty else 0
    counts["is_moonlighting"] = counts["n_groups"] > mode
    return counts.sort_values(["n_groups", "gene"], ascending=[False, True]).reset_index(drop=True), mode


def moonlighting_fraction(members: Iterable[str], moonlighting_ids: set[str]) -> float:
    """Fraction of `members` flagged as moonlighting (present in moonlighting_ids)."""
    members = list(members)
    if not members:
        return np.nan
    return sum(1 for member in members if member in moonlighting_ids) / len(members)
