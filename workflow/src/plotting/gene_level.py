"""
Gene-Level Plotting
===================

Gene-level visualizations for the depletion feature space (DR = depletion rate, DL = depletion lag):
per-group depletion curves, feature-space scatters, and the cluster-merge
review plot used by notebooks/clustering/finalize_gene_clusters.ipynb. Ported
from DIT_HAP_pipeline/workflow/src/subset_visualization.py plus the notebook's
`visualize_cluster_on_feature_space` helper (design doc §7: gene-level plotting
may hardcode genomics assumptions such as the YES0..YES4 LFC columns).

Usage
-----
    from workflow.src.plotting.gene_level import (
        plot_depletion_curves_for_groups,
        plot_groups_on_feature_space,
        visualize_cluster_on_feature_space,
    )
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from typing import Any

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from scipy import stats
from scipy.stats import gaussian_kde

# 4. Local Imports
from workflow.src.plotting.style import AX_HEIGHT, AX_WIDTH

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Raw and fitted LFC value columns for the 5 depletion time points (hardcoded
# genomics assumption — the DIT-HAP assay always has these 5 columns).
RAW_VALUE_COLS = ["YES0", "YES1", "YES2", "YES3", "YES4"]
FITTED_VALUE_COLS = ["YES0_fitted", "YES1_fitted", "YES2_fitted", "YES3_fitted", "YES4_fitted"]

# 10-color palette for final (merged) clusters — assumes <= 10 clusters.
CLUSTER_COLORS = [
    "#dd8369", "#6b99df", "#98a64e", "#64af6d", "#a78bd9",
    "#d57fbd", "#c4954b", "#4bb29c", "#e0788f", "#4aadce",
]


# =============================================================================
# DEPLETION CURVES
# =============================================================================
def plot_depletion_curves_for_given_genes(
    ax: Axes,
    data_df: pd.DataFrame,
    time_points: list[float],
    genes: list[str],
    gene_column: str,
    title: str,
    use_fitted: bool = False,
    **kwargs: Any,
) -> Axes:
    """Plot each gene's depletion curve in gray with the group centroid overlaid in red."""
    value_cols = FITTED_VALUE_COLS if use_fitted else RAW_VALUE_COLS

    for gene in genes:
        values = data_df.query(f"{gene_column} == @gene")[value_cols].values.flatten().tolist()
        ax.plot(time_points, values, color="gray", alpha=0.5, linewidth=0.5, **kwargs)

    subset_df = data_df.query(f"{gene_column} in @genes")
    if subset_df.shape[0] > 0:
        centroid = [subset_df[col].mean() for col in value_cols]
        ax.plot(time_points, centroid, color="red", linewidth=2, label=f"Centroid (n={len(subset_df)})", **kwargs)

    ax.set_title(f"{title}\n(n={len(subset_df)})")
    ax.set_xlabel("Generations")
    ax.set_ylabel("LFC Value")
    ax.grid(True)
    return ax


def plot_depletion_curves_for_groups(
    data_df: pd.DataFrame,
    time_points: list[float],
    group_column: str,
    gene_column: str,
    col_num: int = 4,
    use_fitted: bool = False,
    **kwargs: Any,
) -> Figure:
    """Grid of depletion-curve subplots, one per group value in group_column."""
    groups = sorted(data_df[group_column].unique().tolist())
    row_num = int(np.ceil(len(groups) / col_num))

    fig, axes = plt.subplots(row_num, col_num, figsize=(AX_WIDTH * col_num, AX_HEIGHT * row_num), sharex=True, sharey=True)
    axes = axes.reshape(1, -1) if (row_num == 1 and col_num == 1) else axes.flatten()

    idx = -1
    for idx, group in enumerate(groups):
        group_df = data_df.query(f"{group_column} == @group")
        genes = group_df[gene_column].unique().tolist()
        ax = plot_depletion_curves_for_given_genes(
            ax=axes[idx],
            data_df=data_df,
            time_points=time_points,
            genes=genes,
            gene_column=gene_column,
            title=f"{group_column} {group}",
            use_fitted=use_fitted,
            **kwargs,
        )
        ax.tick_params(axis="both", which="major", labelleft=True, labelbottom=True)

    for j in range(idx + 1, len(axes)):
        fig.delaxes(axes[j])
    return fig


# =============================================================================
# SINGLE-GENE DEPLETION CURVE (DIT-HAP vs gRNA)
# =============================================================================
# Time points (generations) for the two assays, hardcoded in the source notebook
# compare_with_deletion_library.ipynb cell 34. DIT-HAP samples 5 points, gRNA 6.
DIT_HAP_GENERATIONS = [0.0, 2.352, 5.588, 9.104, 12.48]
GRNA_GENERATIONS = [0.0, 4.8, 7.9, 11.4, 14.7, 17.8]
# Raw gRNA per-timepoint LFC columns in HD_gRNA_data.csv (t0..t5).
GRNA_VALUE_COLS = ["t0", "t1", "t2", "t3", "t4", "t5"]


def sigmoid_gompertz(x: np.ndarray, A: float, DR: float, DL: float) -> np.ndarray:
    """Gompertz depletion curve A*exp(-exp(alpha*(DL-x)+1)) with a stable exponent.

    Ported from compare_with_deletion_library.ipynb's sigmoid_function, renamed
    to the release column vocabulary (DR/DL; the notebook used um/lam). A==0
    yields a flat zero curve; the exponent is clipped to [-700, 700] to avoid
    overflow in np.exp.
    """
    if A == 0:
        return np.zeros_like(x)
    alpha = (DR * np.e) / A
    u = alpha * (DL - x) + 1
    exponent = np.clip(u, -700, 700)
    return A * np.exp(-np.exp(exponent))


def plot_gene_depletion_curve(
    ax: Axes,
    dit_row: pd.Series,
    grna_row: pd.Series | None,
    title: str,
    dit_generations: list[float] = DIT_HAP_GENERATIONS,
    grna_generations: list[float] = GRNA_GENERATIONS,
) -> Axes:
    """Plot one gene's DIT-HAP depletion curve (points + Gompertz fit + inflection slope), optionally overlaying gRNA.

    dit_row must carry A/DR/DL plus the RAW_VALUE_COLS (YES0..YES4). grna_row,
    when given, must carry the GRNA_VALUE_COLS (t0..t5); pass None to render the
    DIT-HAP curve alone. title is shown verbatim (the caller passes gene_name).
    """
    A, DR, DL = float(dit_row["A"]), float(dit_row["DR"]), float(dit_row["DL"])
    x = np.linspace(0, 13, 100)
    y_fit = sigmoid_gompertz(x, A, DR, DL)

    # Inflection slope segment: the linear part tangent at the inflection point,
    # spanning the generations where the fitted curve rises (byte-faithful to
    # the notebook's x_slope/y_slope construction).
    xstart = max(DL, 0)
    xend = max(DL + A / DR, 1) if DR != 0 else 1
    x_slope = np.linspace(xstart, xend, 100)
    y_slope = (x_slope - DL) * DR

    y_dit = dit_row[RAW_VALUE_COLS].to_numpy(dtype=float)

    ax.plot(x, y_fit, label="Fitted Curve", linestyle="--", alpha=0.5)
    ax.plot(x_slope, y_slope, label="Inflection Point", linestyle=":", alpha=0.5)
    ax.plot(dit_generations, y_dit, marker="o", label=f"DIT_HAP\n(DR={DR:.3g})")

    if grna_row is not None:
        y_grna = grna_row[GRNA_VALUE_COLS].to_numpy(dtype=float)
        ax.plot(grna_generations, y_grna, marker="o", label="gRNA")

    ax.set_ylim(-2, 9)
    ax.set_title(title)
    ax.set_xlabel("Generation")
    ax.set_ylabel("LFC")
    ax.legend()
    return ax


# =============================================================================
# FEATURE SPACE
# =============================================================================
def plot_given_genes_on_feature_space(
    ax: Axes,
    data_df: pd.DataFrame,
    genes: list[str],
    gene_column: str,
    title: str,
    x_feature: str = "DR",
    y_feature: str = "DL",
    cmap: str | LinearSegmentedColormap = "viridis",
    label: str = "Selected Genes",
    title_with_count: bool = True,
    **kwargs: Any,
) -> Axes:
    """Scatter all genes in light gray, highlighting the given subset in color."""
    ax.scatter(data_df[x_feature], data_df[y_feature], color="lightgray", alpha=0.4, **kwargs)

    subset_df = data_df.query(f"`{gene_column}` in @genes")
    x_subset = subset_df[x_feature]
    y_subset = subset_df[y_feature]
    try:
        xy_subset = np.vstack([x_subset, y_subset])
        gaussian_kde(xy_subset)(xy_subset)  # density check preserved from source
        ax.scatter(x_subset, y_subset, color=cmap, **kwargs, label=f"{label} (n={len(subset_df)})")
    except Exception:
        ax.scatter(x_subset, y_subset, color="red", **kwargs, label=f"{label} (n={len(subset_df)})")

    ax.set_title(f"{title}\n(n={len(subset_df)})" if title_with_count else f"{title}")
    ax.set_xlabel(x_feature)
    ax.set_ylabel(y_feature)
    return ax


def plot_groups_on_feature_space(
    data_df: pd.DataFrame,
    group_column: str,
    gene_column: str,
    col_num: int = 4,
    x_feature: str = "DR",
    y_feature: str = "DL",
    **kwargs: Any,
) -> Figure:
    """Grid of feature-space subplots, one per group value in group_column."""
    groups = sorted(data_df[group_column].unique().tolist())
    row_num = int(np.ceil(len(groups) / col_num))

    fig, axes = plt.subplots(row_num, col_num, figsize=(AX_WIDTH * col_num, AX_HEIGHT * row_num), sharex=True, sharey=True)
    axes = axes.reshape(1, -1) if (row_num == 1 and col_num == 1) else axes.flatten()

    idx = -1
    for idx, group in enumerate(groups):
        group_df = data_df.query(f"{group_column} == @group")
        genes = group_df[gene_column].unique().tolist()
        ax = plot_given_genes_on_feature_space(
            ax=axes[idx],
            data_df=data_df,
            genes=genes,
            gene_column=gene_column,
            title=f"{group_column} {group}",
            x_feature=x_feature,
            y_feature=y_feature,
            **kwargs,
        )
        ax.tick_params(axis="both", which="major", labelleft=True, labelbottom=True)

    for j in range(idx + 1, len(axes)):
        fig.delaxes(axes[j])
    return fig


# =============================================================================
# CLUSTER-MERGE REVIEW (used by finalize_gene_clusters notebook)
# =============================================================================
# 64-color palette sized to match the candidate n_clusters=64 scatter.
MULTI_COLORS = [
    "#452062", "#8be93d", "#6d2bde", "#5bdb62", "#cc3ce7", "#cfe044", "#5447d1", "#e5c83b",
    "#9643bf", "#6cac2e", "#d93eb4", "#63e3a1", "#49288c", "#b8e181", "#6474df", "#da9b33",
    "#c875d8", "#55a95c", "#dd3d86", "#69decb", "#e63923", "#63c0df", "#db632f", "#6897da",
    "#c1712e", "#425594", "#9ca13e", "#913379", "#477126", "#de4867", "#3c8b67", "#cc403c",
    "#afddd8", "#271c45", "#e4cc7f", "#9379c4", "#90782c", "#d9a5df", "#2b4826", "#dc76ab",
    "#d0e2aa", "#321521", "#e1d2bf", "#182f3b", "#d49968", "#537d9a", "#8a3e1f", "#b6bbdc",
    "#652524", "#97bc8f", "#9e2e4c", "#649f9e", "#612643", "#788a63", "#b15f6a", "#385f5d",
    "#e28c7e", "#3f311f", "#d4a4b1", "#605320", "#8f6c8b", "#b09a7e", "#52485f", "#8a6653",
]


def create_gradient_colormap(color: str, name: str) -> LinearSegmentedColormap:
    """Build a white-to-`color` gradient colormap for density-shaded scatter."""
    return LinearSegmentedColormap.from_list(name, ["white", color], N=256)


def plot_cluster_on_axis(
    ax: Axes,
    df: pd.DataFrame,
    cluster_col: str,
    cluster_minus_one: bool = False,
    show_box: bool = True,
    title: str | None = None,
    point_size: int = 8,
    label_fontsize: int = 14,
) -> Axes:
    """Plot one clustering's assignments on a single axis (flat per-cluster color).

    The single-axis analogue of `visualize_cluster_on_feature_space`'s flat-color
    panel, so a multi-variant grid reuses the exact same color mapping (CLUSTER_COLORS,
    modulo-indexed) and DR/DL framing. `cluster_minus_one=True` when cluster ids are
    1-based (final 1..9) but the color list is 0-indexed.
    """
    for cluster, cluster_df in df.groupby(cluster_col, sort=True):
        cluster = int(cluster)
        x = cluster_df["DR"]
        y = cluster_df["DL"]
        cluster_idx = cluster - 1 if cluster_minus_one else cluster
        cluster_color = CLUSTER_COLORS[cluster_idx % len(CLUSTER_COLORS)]
        ax.scatter(x, y, c=cluster_color, s=point_size, label=f"Cluster {cluster} (n={len(x)})")
        centroid_x, centroid_y = x.mean(), y.mean()
        if show_box:
            ax.text(centroid_x, centroid_y, f"{cluster}", fontweight="bold", fontsize=label_fontsize,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.8))
        else:
            ax.text(centroid_x, centroid_y, f"{cluster}", fontweight="bold", fontsize=label_fontsize)

    ax.set_xlabel("DR")
    ax.set_ylabel("DL")
    ax.set_xlim(-0.15, 1.45)
    if title:
        ax.set_title(title)
    return ax


def visualize_cluster_on_feature_space(
    df: pd.DataFrame,
    cluster_col: str,
    cluster_minus_one: bool = False,
    legend: bool = False,
    show_box: bool = False,
) -> Figure:
    """Two-panel feature-space scatter (KDE-shaded + flat color) with per-cluster labels.

    `cluster_minus_one=True` when cluster ids are 1-based (cluster is 1..9)
    but the color lists are 0-indexed.
    """
    fig, axes = plt.subplots(1, 2, figsize=(18, 7), sharex=True, sharey=True)
    for cluster, cluster_df in df.groupby(cluster_col, sort=True):
        cluster = int(cluster)
        x = cluster_df["DR"]
        y = cluster_df["DL"]
        cluster_idx = cluster - 1 if cluster_minus_one else cluster

        # Index palettes with modulo so the 64-cluster candidate review does not
        # overrun the 10-color merged palette. For the merged 1..9 case (idx 0..8)
        # this is a no-op, so final-cluster colors stay byte-faithful.
        multi_color = MULTI_COLORS[cluster_idx % len(MULTI_COLORS)]
        cluster_color = CLUSTER_COLORS[cluster_idx % len(CLUSTER_COLORS)]

        try:
            xy = np.vstack([x, y])
            z = stats.gaussian_kde(xy)(xy)
            cmap = create_gradient_colormap(multi_color, f"{cluster}_gradient")
            axes[0].scatter(x, y, c=z, s=20, cmap=cmap, label=f"Cluster {cluster} (n={len(x)})")
        except Exception:
            axes[0].scatter(x, y, c=multi_color, s=20, label=f"Cluster {cluster}\n(n={len(x)})")

        axes[1].scatter(x, y, c=cluster_color, s=20, label=f"Cluster {cluster}\n(n={len(x)})")
        centroid_x, centroid_y = x.mean(), y.mean()
        for ax in axes:
            if show_box:
                ax.text(centroid_x, centroid_y, f"{cluster}", fontweight="bold", fontsize=22,
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.8))
            else:
                ax.text(centroid_x, centroid_y, f"{cluster}", fontweight="bold", fontsize=22)

    for ax in axes:
        ax.set_xlabel("DR")
        ax.set_ylabel("DL")
        ax.tick_params(axis="both", which="major", labelsize=20, labelleft=True, labelbottom=True)
        if legend:
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=16, frameon=True, markerscale=5)
        ax.set_xlim(-0.15, 1.45)

    plt.tight_layout()
    return fig
