"""Generic histogram grid rendering.

Supersedes ``distribution.py`` (one panel per metric, raw values) and
``read_counts.py`` (pre-binned counts replayed via weights, with a cutoff marker
and a figure-level footer). Both are grids of histograms, but their inputs
differ in form, so each has its own entry point here.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.1.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence
from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from loguru import logger
from matplotlib.axes import Axes

from workflow.src.figures import (
    FURNITURE_COLOR,
    PanelShape,
    apply_house_style,
    apply_log_scale,
    fit_panels,
    grid_axes,
    panel_labels,
    save_dual,
)

from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================
# cns.histplot addresses data by column name, so an incoming Series is boxed
# under this name. It never reaches the axis labels, which callers supply.
_VALUE_COLUMN = "value"


# =============================================================================
# CORE LOGIC — single-axes primitive
# =============================================================================
def draw_histogram_panel(
    ax: Axes,
    data: pd.Series,
    *,
    bins: int,
    log_scale: bool = False,
    xlabel: str = "Value",
    ylabel: str = "Frequency",
    title: str = "",
    bin_edges: np.ndarray | None = None,
) -> None:
    """Histogram one Series onto ax; the shared primitive behind both grid renderers.

    Handles positive-value filtering in log mode, the empty-data placeholder,
    and the house bar style. When ``bin_edges`` is given it is used as-is (for
    callers coordinating ranges across panels); otherwise edges come from the
    data itself, laid evenly in log10 space under ``log_scale``.
    """
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    values = data.dropna()
    if log_scale:
        values = values[values > 0]

    if values.empty:
        logger.warning("Panel has no valid data")
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
        return

    if bin_edges is not None:
        edges = bin_edges
    elif log_scale:
        edges = np.logspace(np.log10(values.min()), np.log10(values.max()), bins + 1)
    else:
        edges = bins

    # cns.histplot takes a frame and a column name, and the Series may be
    # unnamed or share a name with nothing here, so it is boxed under a fixed
    # internal column. The scale is set afterwards rather than through
    # log_scale=: the edges are already laid out in log10 space above, and
    # passing both would have seaborn re-bin them.
    cns.histplot(data=values.to_frame(_VALUE_COLUMN), x=_VALUE_COLUMN, bins=edges, ax=ax, shrink=0.9)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if log_scale:
        apply_log_scale(ax, x=True, y=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_histogram_grid_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    value_columns: Sequence[str],
    bins: int,
    xlabel: str = "Value",
    ylabel: str = "Frequency",
    show_summary_stats: bool = True,
    n_cols: int = 4,
    share_y_range: bool = False,
    shape: PanelShape = PanelShape.SQUARE,
) -> None:
    """Render one histogram panel per named value column in a fixed-pitch grid."""
    logger.info("Rendering histogram grid figure...")

    require_columns(df, value_columns, context="histogram input")

    if df.empty or not value_columns:
        logger.warning("No data to plot!")
        return
    
    if n_cols < 1:
        raise ValueError(f"n_cols must be at least 1, got {n_cols}")

    apply_house_style()

    n_panels = len(value_columns)
    n_rows = (n_panels + n_cols - 1) // n_cols
    logger.info(f"Creating figure with {n_rows}x{n_cols} grid ({n_panels} panels)...")

    labels = panel_labels(n_panels)
    axes = grid_axes(n_rows, n_cols, labels=labels, shape=shape)

    finite_maxes = [m for c in value_columns if (m := df[c].max()) == m]
    shared_y_max = max(finite_maxes, default=0.0) if share_y_range else None

    for ax in axes[n_panels:]:
        ax.set_visible(False)

    for index, (label, column) in enumerate(zip(labels, value_columns, strict=True)):
        ax = axes[index]
        data = df[column].dropna()
        logger.info(f"  Panel {label}: {column} (n={len(data)})")

        had_data = data.notna().any()
        draw_histogram_panel(ax, data, bins=bins, xlabel=xlabel, ylabel=ylabel, title=column)
        if had_data and show_summary_stats:
            stats_text = f"n = {len(data):,}\nMean = {data.mean():.3f}\nStd = {data.std():.3f}"
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, verticalalignment="top")
        if shared_y_max is not None and had_data:
            ax.set_ylim(0, shared_y_max)

    fit_panels()

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")


@logger.catch(reraise=True)
def render_grouped_histogram_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    value_column: str,
    row_key: str,
    col_key: str,
    bins: int = 50,
    log_scale: bool = False,
    xlabel: str = "Value",
    ylabel: str = "Frequency",
    marker_value: float | None = None,
    marker_label: str = "",
    marker_on_col_value: str | None = None,
    upper_quantile: float | None = None,
    share_x_range: bool = True,
    share_y_range: bool = False,
    shape: PanelShape = PanelShape.SQUARE,
) -> None:
    """Histogram one value column into one panel per row_key/col_key group in a fixed-pitch grid.

    With ``log_scale`` the values must be strictly positive: bin edges are laid
    evenly in log10 space via ``np.logspace`` and the axis becomes a true log
    scale, so bars are geometrically identical to pre-transformed linear binning
    while ticks show real read counts. The cutoff marker is then placed at the
    raw value (e.g. x=8), not at its log.

    ``upper_quantile`` (e.g. 0.999999) drops values above each group's own
    quantile before binning — QC read-count tables carry one or two astronomic
    outliers that otherwise stretch the log axis across empty decades.

    ``share_x_range`` puts every panel on common bin edges so bars are
    comparable left-to-right; ``share_y_range`` is off by default because group
    frequencies differ by orders of magnitude and one global top flattens the
    smaller panels.

    The row_key value rides in the first column's ylabel rather than the title.
    """
    require_columns(df, [row_key, col_key, value_column], context="grouped histogram input")

    if df.empty:
        logger.warning("No values to plot!")
        return

    grouped = df.groupby([row_key, col_key], sort=True)
    if grouped.ngroups == 0:
        logger.warning("No groups to plot!")
        return

    def _panel_values(group_df: pd.DataFrame) -> pd.Series:
        data = group_df[value_column].dropna()
        if log_scale:
            data = data[data > 0]
        if upper_quantile is not None and not data.empty:
            threshold = data.quantile(upper_quantile)
            dropped = int((data > threshold).sum())
            if dropped:
                logger.info(f"upper_quantile={upper_quantile}: dropping {dropped} outlier(s) above {threshold:.4g}")
            data = data[data <= threshold]
        return data

    panel_values = {key: _panel_values(group_df) for key, group_df in grouped}
    if all(data.empty for data in panel_values.values()):
        logger.warning("All panels are empty after filtering!")
        return

    apply_house_style()

    n_cols = df[col_key].nunique()
    n_rows = df[row_key].nunique()

    row_values = sorted(df[row_key].unique())
    col_values = sorted(df[col_key].unique())

    labels = panel_labels(n_rows * n_cols)
    axes = grid_axes(n_rows, n_cols, labels=labels, shape=shape)

    # A row/col combination absent from the data keeps its cell but shows no
    # frame, so the remaining panels stay in their own row and column.
    for r, row_val in enumerate(row_values):
        for c, col_val in enumerate(col_values):
            if (row_val, col_val) not in panel_values:
                axes[r * n_cols + c].set_visible(False)

    if share_x_range and log_scale:
        lo = min(d.min() for d in panel_values.values() if not d.empty)
        hi = max(d.max() for d in panel_values.values() if not d.empty)
        shared_edges = np.logspace(np.log10(lo), np.log10(hi), bins + 1)
    else:
        shared_edges = None

    for (row_val, col_val), data in panel_values.items():
        r = row_values.index(row_val)
        c = col_values.index(col_val)
        ax = axes[r * n_cols + c]

        panel_ylabel = f"{row_val}\n{ylabel}" if c == 0 else ylabel
        draw_histogram_panel(
            ax, data,
            bins=bins,
            log_scale=log_scale,
            xlabel=xlabel,
            ylabel=panel_ylabel,
            title=str(col_val),
            bin_edges=shared_edges,
        )

        if share_y_range and not data.empty:
            y_max = max(
                int(np.histogram(d.to_numpy(), bins=(shared_edges if shared_edges is not None else bins))[0].max())
                for d in panel_values.values()
                if not d.empty
            )
            ax.set_ylim(0, y_max * 1.05)

        logger.info(f"  Panel ({r},{c}): {row_val} {col_val} (n={len(data):,})")

        if marker_value is None:
            continue
        if marker_on_col_value is not None and col_val != marker_on_col_value:
            continue

        ax.axvline(marker_value, color=FURNITURE_COLOR, linestyle="--", label=marker_label)
        ax.legend(loc="upper right")

    fit_panels()

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
