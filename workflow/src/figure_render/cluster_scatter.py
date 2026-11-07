"""
Cluster Scatter Rendering
=========================

DR/DL cluster-assignment figures in the house cnsplots style: the per-variant
review page (density + flat colour) and the all-variants grid.

The legacy ``plotting.gene_level.visualize_cluster_on_feature_space`` remains
for ``notebooks/clustering/finalize_gene_clusters.ipynb``, whose headless run
uses the computing environment where cnsplots is not installed.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-07
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Mapping, Sequence

from loguru import logger
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from figures import (
    CLUSTER_GRID_SCATTER_S,
    PanelShape,
    apply_house_style,
    cluster_color,
    fit_panels,
    grid_axes,
    panel_labels,
)
from plotting.gene_level import DR_XLIM, create_gradient_colormap

# =============================================================================
# CONSTANTS
# =============================================================================
# Right-hand strip reserved for the shared cluster legend, as a tight_layout
# rect and the matching axes-fraction anchor (same split the scatter renderers use).
_LEGEND_STRIP_LEFT = 0.86
LEGEND_ANCHOR = (_LEGEND_STRIP_LEFT, 0.5)
LEGEND_LAYOUT_RECT = (0, 0, _LEGEND_STRIP_LEFT, 1)

# Space kept above the panels for a suptitle.
_TITLE_TOP = 0.95


def _layout_rect(*, legend: bool, title: bool) -> tuple[float, float, float, float] | None:
    """The tight_layout rect reserving the legend strip and/or the suptitle band."""
    if not legend and not title:
        return None
    return (0, 0, _LEGEND_STRIP_LEFT if legend else 1, _TITLE_TOP if title else 1)


# =============================================================================
# CORE LOGIC
# =============================================================================
def _cluster_colors(cluster_ids: Sequence[int], *, one_based: bool) -> dict[int, str]:
    """House colour per cluster id; 0-based ids start at the first palette colour."""
    offset = 0 if one_based else 1
    return {int(cluster): cluster_color(int(cluster) + offset) for cluster in cluster_ids}


def _chip(ax: Axes, cluster: int, x: float, y: float, *, show_box: bool) -> None:
    """Write the cluster id at its centroid, optionally on a white chip."""
    bbox = (
        {
            "boxstyle": "round,pad=0.15",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.8,
        }
        if show_box
        else None
    )
    ax.annotate(f"{cluster}", (x, y), fontweight="bold", ha="center", va="center", bbox=bbox)


def _cluster_ids(df: pd.DataFrame) -> list[int]:
    """Sorted integer cluster ids, so palette assignment and drawing order agree."""
    return sorted(int(cluster) for cluster in df["cluster"].unique())


def scatter_clusters(
    ax: Axes, df: pd.DataFrame, colors: Mapping[int, str], *, show_box: bool = True
) -> None:
    """Draw one DR/DL scatter on one axes: per-cluster colour plus a centroid id chip."""
    for cluster in _cluster_ids(df):
        cluster_df = df.loc[df["cluster"] == cluster]
        ax.scatter(
            cluster_df["DR"],
            cluster_df["DL"],
            c=colors[cluster],
            s=CLUSTER_GRID_SCATTER_S,
            rasterized=True,
        )
        _chip(
            ax,
            cluster,
            cluster_df["DR"].mean(),
            cluster_df["DL"].mean(),
            show_box=show_box,
        )
    ax.set_xlabel("DR")
    ax.set_ylabel("DL")


def _draw_density_panel(
    ax: Axes, frame: pd.DataFrame, colors: Mapping[int, str], *, show_box: bool
) -> None:
    """Colour each cluster's points by its own 2D density, falling back to flat colour."""
    for cluster in _cluster_ids(frame):
        cluster_df = frame.loc[frame["cluster"] == cluster]
        x = cluster_df["DR"]
        y = cluster_df["DL"]
        try:
            density = stats.gaussian_kde(np.vstack([x, y]))(np.vstack([x, y]))
            cmap = create_gradient_colormap(colors[cluster], f"{cluster}_gradient")
            ax.scatter(x, y, c=density, s=CLUSTER_GRID_SCATTER_S, cmap=cmap, rasterized=True)
        except (np.linalg.LinAlgError, ValueError):
            ax.scatter(x, y, c=colors[cluster], s=CLUSTER_GRID_SCATTER_S, rasterized=True)
        _chip(ax, cluster, x.mean(), y.mean(), show_box=show_box)
    ax.set_xlabel("DR")
    ax.set_ylabel("DL")
    ax.set_xlim(DR_XLIM)


def _add_cluster_legend(fig: Figure, frame: pd.DataFrame, colors: Mapping[int, str]) -> None:
    """One figure-level legend entry per cluster, with its gene count."""
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            markersize=4,
            markerfacecolor=colors[int(cluster)],
            markeredgecolor="none",
            label=f"Cluster {int(cluster)} (n={int((frame['cluster'] == cluster).sum())})",
        )
        for cluster in sorted(colors)
    ]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=LEGEND_ANCHOR, frameon=False)


def build_cluster_scatter_figure(
    df: pd.DataFrame,
    *,
    cluster_col: str,
    one_based: bool,
    show_box: bool = True,
    legend: bool = False,
    title: str | None = None,
) -> Figure:
    """Build the two-panel (density + flat colour) DR/DL page for one clustering."""
    logger.info(f"Rendering cluster scatter for {len(df)} genes...")
    apply_house_style()

    # The clustered frames also carry a `cluster` column when the page is the
    # pre-merge raw_cluster one, so overwrite rather than rename into a duplicate.
    frame = df[["DR", "DL", cluster_col]].rename(columns={cluster_col: "cluster"})
    frame["cluster"] = frame["cluster"].astype(int)
    colors = _cluster_colors(frame["cluster"].unique(), one_based=one_based)

    axes = grid_axes(1, 2, shape=PanelShape.SQUARE, share_x=True, share_y=True)
    _draw_density_panel(axes[0], frame, colors, show_box=show_box)
    scatter_clusters(axes[1], frame, colors, show_box=show_box)

    fig = plt.gcf()
    if legend:
        _add_cluster_legend(fig, frame, colors)
    if title is not None:
        fig.suptitle(title)
    fit_panels(rect=_layout_rect(legend=legend, title=title is not None))
    return fig


def build_cluster_grid_figure(
    frames: Mapping[str, pd.DataFrame],
    *,
    col_num: int = 3,
    title: str | None = None,
) -> Figure:
    """Build one aligned DR/DL panel per variant, sharing cluster colours across panels."""
    labels = list(frames)
    logger.info(f"Rendering cluster grid with {len(labels)} variants...")
    apply_house_style()

    n_cols = min(col_num, len(labels))
    n_rows = -(-len(labels) // n_cols)
    axes = grid_axes(n_rows, n_cols, labels=panel_labels(len(labels)), share_x=True, share_y=True)

    normalized = {
        label: frame[["DR", "DL", "cluster"]].astype({"cluster": int})
        for label, frame in frames.items()
    }
    cluster_ids = sorted({int(c) for frame in normalized.values() for c in frame["cluster"]})
    colors = _cluster_colors(cluster_ids, one_based=True)

    for ax, (label, frame) in zip(axes, normalized.items(), strict=False):
        scatter_clusters(ax, frame, colors)
        ax.set_title(label)
    for ax in axes[len(labels) :]:
        ax.set_visible(False)

    fig = plt.gcf()
    if title is not None:
        fig.suptitle(title)
    fit_panels(rect=_layout_rect(legend=False, title=title is not None))
    return fig
