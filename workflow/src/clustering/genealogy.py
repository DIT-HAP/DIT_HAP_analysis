"""
Cross-k cluster genealogy
=========================

How one k's clusters become the next k's. Sweeping k answers "how many", but it
does not say WHERE the extra clusters come from: a cluster count whose members
are all fragments of one k=5 cluster is buying refinements of a tail, not new
groups.

The genealogy is computed by overlap, not by any tree structure: each cluster at
k+1 is matched to the cluster at k that holds most of its genes. Only the
dominant parent is kept, because a cluster whose genes are split across several
parents would turn the alluvial into a lattice. A cluster with no parent meeting
the minimum share starts its own lineage.

Input
-----
- grid_labels.parquet (method, k, gene, cluster), the DR-ascending cluster ids.

Output
------
- A long table, one row per parent->child link above the minimum overlap.

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-10-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from collections.abc import Sequence

# 2. Data Processing Imports
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle

# 4. Local Imports
from plotting.gene_level import CLUSTER_COLORS


def cluster_color(cluster: int) -> str:
    """One colour per cluster ID within its column (same mapping as the scatter pages)."""
    return CLUSTER_COLORS[(int(cluster) - 1) % len(CLUSTER_COLORS)]


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# One row per parent->child link; `k` is the PARENT's, so a row reads "cluster a
# at k became cluster b at k+1".
GENEALOGY_COLUMNS: tuple[str, ...] = (
    "k",
    "cluster",
    "n_genes",
    "child_cluster",
    "child_genes",
    "overlap",
    "overlap_share_of_child",
    "n_parents",
)

# Below this share of the child's genes, a parent link is not drawn: a handful
# of genes otherwise draw a band across the whole figure (v2 MIN_LINK_SHARE).
MIN_LINK_SHARE = 0.05

# Alluvial geometry (v2 render_genealogy_panel values): bar width in k units,
# the column share a bar must reach before it is worth numbering, its type size,
# and the band alpha (barely-filled bands so 15 columns of them don't pile up).
BAR_HALF_WIDTH = 0.16
LABEL_MIN_SHARE = 0.07
LABEL_FONTSIZE = 5
BAND_ALPHA = 0.25


# =============================================================================
# CORE LOGIC — table
# =============================================================================
def cluster_genealogy(labels: pd.DataFrame, k_values: Sequence[int]) -> pd.DataFrame:
    """Link every k's clusters of ONE method to the previous k's by gene overlap."""
    # The gene universe, fixed once: parent and child arrays must sit in the
    # SAME gene order, and each gene exactly once (a method frame carries all
    # its k's rows, so its gene column repeats per k).
    universe = sorted(labels["gene"].unique())
    partitions: dict[int, np.ndarray] = {}
    for k in k_values:
        subset = labels[labels["k"] == k].set_index("gene")
        missing = set(universe) - set(subset.index)
        if missing:
            raise ValueError(f"Genealogy: {len(missing)} genes missing at k={k}")
        partitions[k] = subset.loc[universe, "cluster"].to_numpy()

    rows = []
    for parent_k, child_k in zip(k_values, k_values[1:], strict=False):
        parent = partitions[parent_k]
        child = partitions[child_k]
        # One row per (parent, child) pair that actually shares genes.
        pairs = pd.crosstab(pd.Series(parent, name="cluster"), pd.Series(child, name="child_cluster"))
        child_sizes = pairs.sum(axis=0)
        for child_label in pairs.columns:
            column = pairs[child_label]
            shared = column[column > 0]
            share = float(shared.max() / child_sizes[child_label])
            if share < MIN_LINK_SHARE:
                continue
            rows.append(
                {
                    "k": parent_k,
                    "cluster": int(shared.idxmax()),
                    "n_genes": int((parent == shared.idxmax()).sum()),
                    "child_cluster": int(child_label),
                    "child_genes": int(child_sizes[child_label]),
                    "overlap": int(shared.max()),
                    "overlap_share_of_child": share,
                    "n_parents": int(len(shared)),
                }
            )
    return pd.DataFrame(rows).reindex(columns=GENEALOGY_COLUMNS)


# =============================================================================
# CORE LOGIC — alluvial panel
# =============================================================================
def render_genealogy_panel(ax: Axes, genealogy: pd.DataFrame, *, k_values: Sequence[int], method: str) -> None:
    """Render how one method's clusters split from k to k, as a stacked alluvial."""
    # Each column is one k; a column's bars are that partition's clusters, sized
    # by their share of the genes, coloured by the cluster they descend from at
    # the smallest k. Reading it across, a bar that stays one band survived; a
    # bar that fans into several is the extra k buying a refinement of it.
    if genealogy.empty:
        ax.set(xlabel="k", ylabel="Share of genes", title=f"Genealogy ({method})", xlim=(-0.5, len(k_values) - 0.5), ylim=(0, 1))
        ax.text(0.5, 0.5, "no stable links", transform=ax.transAxes, ha="center", va="center")
        return

    # Cluster ids are only unique WITHIN a k, so every lookup is keyed by
    # (k, cluster). Bars and bands take the CLUSTER's own id colour (the same
    # per-id mapping as the scatter pages), so at one k each cluster reads in
    # its own hue and a split visibly hands over to the children's colours.
    x_lookup = {k: position for position, k in enumerate(k_values)}

    # Columns are ordered by where each cluster's parent sat in the previous
    # column, NOT by cluster id — id order makes sibling bands cross the column.
    members: dict[int, set[int]] = {k: set() for k in k_values}
    parent_of: dict[tuple[int, int], int] = {}
    for _, row in genealogy.iterrows():
        parent_k = int(row["k"])
        members[parent_k].add(int(row["cluster"]))
        members[parent_k + 1].add(int(row["child_cluster"]))
        # First row wins per child (the table's dominant-parent rule).
        parent_of.setdefault((parent_k + 1, int(row["child_cluster"])), int(row["cluster"]))

    order: dict[int, list[int]] = {}
    previous_rank: dict[int, int] = {}
    for k in k_values:
        ranked = sorted(members[k], key=lambda c: (previous_rank.get(parent_of.get((k, c), -1), 0), c))
        order[k] = ranked
        previous_rank = {c: rank for rank, c in enumerate(ranked)}

    bars: dict[tuple[int, int], tuple[float, float]] = {}
    for k in k_values:
        cursor = 0.0
        for cluster in order[k]:
            own = genealogy[(genealogy["k"] == k) & (genealogy["cluster"] == cluster)]
            size = (
                float(own["n_genes"].iloc[0]) if not own.empty
                else float(genealogy[(genealogy["k"] == k - 1) & (genealogy["child_cluster"] == cluster)]["child_genes"].iloc[0])
            )
            bars[(k, cluster)] = (cursor, cursor + size)
            cursor += size

    # Band stubs inside each bar, stacked in the order the bars are drawn — any
    # other order makes two bands cross inside the bar they share.
    total = max(high for _, high in bars.values())

    def scale(value: float) -> float:
        """A y position in gene counts -> axis units."""
        return value / total

    for parent_k, child_k in zip(k_values, k_values[1:], strict=False):
        links = genealogy[genealogy["k"] == parent_k]
        rank = {c: position for position, c in enumerate(order[child_k])}
        outgoing_cursor = {key: low for key, (low, _) in bars.items() if key[0] == parent_k}
        incoming_cursor = {key: low for key, (low, _) in bars.items() if key[0] == child_k}
        for parent in sorted(order[parent_k], key=lambda c: previous_rank.get(c, 0)):
            children = (
                links[links["cluster"] == parent]
                .sort_values(by="child_cluster", key=lambda column: column.map(rank))
            )
            for _, row in children.iterrows():
                child = int(row["child_cluster"])
                start = outgoing_cursor[(parent_k, parent)]
                child_start = incoming_cursor[(child_k, child)]
                steps = np.linspace(x_lookup[parent_k], x_lookup[child_k], 16)
                # Cosine easing: bands leave and arrive vertically, so a split
                # reads as a fork rather than a diagonal crossing.
                weight = (1 - np.cos(np.pi * (steps - steps[0]) / (steps[-1] - steps[0]))) / 2
                ax.fill_between(
                    steps,
                    scale(start + (child_start - start) * weight),
                    scale(start + float(row["child_genes"]) + (child_start + float(row["n_genes"]) - start - float(row["child_genes"])) * weight),
                    color=cluster_color(child),
                    alpha=BAND_ALPHA, linewidth=0,
                )
                outgoing_cursor[(parent_k, parent)] = start + float(row["child_genes"])
                incoming_cursor[(child_k, child)] = child_start + float(row["n_genes"])

    for (k, cluster), (low, high) in bars.items():
        ax.add_patch(
            Rectangle(
                (x_lookup[k] - BAR_HALF_WIDTH, scale(low)), 2 * BAR_HALF_WIDTH, scale(high - low),
                facecolor=cluster_color(cluster),
                edgecolor="none",
            )
        )
        # Only bars with room for a label get one; at k=15 most are slivers the
        # reader takes by position.
        if high - low > LABEL_MIN_SHARE * total:
            ax.text(x_lookup[k], scale((low + high) / 2), str(cluster), ha="center", va="center", fontsize=LABEL_FONTSIZE)

    ax.set(
        xlabel="k", ylabel="Share of genes",
        title=f"Genealogy ({method})",
        xticks=list(x_lookup.values()), xticklabels=[str(k) for k in x_lookup],
        xlim=(-0.5, len(k_values) - 0.5), ylim=(0, 1.0),
    )
    ax.invert_yaxis()  # the largest cluster at the top, where the eye starts