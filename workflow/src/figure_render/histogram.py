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
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
import sys
from typing import Any

# Add workflow/src to path for imports
_SCRIPT_DIR = Path(__file__).resolve().parent
_SRC_DIR = _SCRIPT_DIR.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

import cnsplots as cns
from loguru import logger
from matplotlib.axes import Axes
import numpy as np
import pandas as pd

from figures import (
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

# Stand-in column for the un-faceted call (col_key=None). Carries one blank
# value, so the grid collapses to a single column and the row key alone labels
# each panel.
_SINGLE_COLUMN = "_single"


# =============================================================================
# CORE LOGIC — single-axes primitive
# =============================================================================
def draw_histogram_panel(
    ax: Axes,
    data: pd.Series,
    *,
    bins: int | np.ndarray,
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
    cns.histplot(
        data=values.to_frame(_VALUE_COLUMN), x=_VALUE_COLUMN, bins=edges, ax=ax, shrink=0.9
    )
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


def _panel_values(
    group_df: pd.DataFrame,
    *,
    value_column: str,
    log_scale: bool,
    upper_quantile: float | None,
) -> pd.Series:
    """Drop nulls, non-positive values under log scale, and per-group upper outliers."""
    data = group_df[value_column].dropna()
    if log_scale:
        data = data[data > 0]
    if upper_quantile is None or data.empty:
        return data

    threshold = data.quantile(upper_quantile)
    dropped = int((data > threshold).sum())
    if dropped:
        logger.info(
            f"upper_quantile={upper_quantile}: dropping {dropped} outlier(s) above {threshold:.4g}"
        )
    return data[data <= threshold]


def _collect_panel_values(
    groups: Iterable[tuple[Any, pd.DataFrame]],
    *,
    value_column: str,
    log_scale: bool,
    upper_quantile: float | None,
) -> dict[Any, pd.Series]:
    """Filter every row_key/col_key group down to the values its panel will bin."""
    return {
        key: _panel_values(
            group_df, value_column=value_column, log_scale=log_scale, upper_quantile=upper_quantile
        )
        for key, group_df in groups
    }


def _hide_absent_cells(
    axes: Sequence[Axes],
    panel_values: Mapping[Any, pd.Series],
    row_values: Sequence[Any],
    col_values: Sequence[Any],
    n_cols: int,
) -> None:
    """Hide the frame of every row/col cell absent from the data, keeping its slot."""
    for r, row_val in enumerate(row_values):
        for c, col_val in enumerate(col_values):
            if (row_val, col_val) not in panel_values:
                axes[r * n_cols + c].set_visible(False)


def _shared_log_edges(
    panel_values: Mapping[Any, pd.Series],
    *,
    share_x_range: bool,
    log_scale: bool,
    bins: int | np.ndarray,
) -> np.ndarray | None:
    """Common log-spaced bin edges across panels, or None when ranges stay per-panel."""
    if not (share_x_range and log_scale):
        return None

    non_empty = [data for data in panel_values.values() if not data.empty]
    lo = min(data.min() for data in non_empty)
    hi = max(data.max() for data in non_empty)
    return np.logspace(np.log10(lo), np.log10(hi), bins + 1)


def _shared_y_top(
    panel_values: Mapping[Any, pd.Series],
    *,
    share_y_range: bool,
    shared_edges: np.ndarray | None,
    bins: int | np.ndarray,
) -> float | None:
    """Highest bar count across panels, with headroom, or None without a shared y range."""
    if not share_y_range:
        return None

    edges = shared_edges if shared_edges is not None else bins
    counts = [
        int(np.histogram(data.to_numpy(), bins=edges)[0].max())
        for data in panel_values.values()
        if not data.empty
    ]
    return max(counts) * 1.05


def _grouped_panel_titles(
    row_val: Any, col_val: Any, c: int, *, un_faceted: bool, ylabel: str
) -> tuple[str, str]:
    """Return the panel ylabel/title, carrying the row value in the first column."""
    if un_faceted:
        return ylabel, str(row_val)

    panel_ylabel = f"{row_val}\n{ylabel}" if c == 0 else ylabel
    return panel_ylabel, str(col_val)


def _draw_grouped_panel(
    ax: Axes,
    data: pd.Series,
    *,
    panel_ylabel: str,
    panel_title: str,
    bins: int | np.ndarray,
    log_scale: bool,
    xlabel: str,
    bin_edges: np.ndarray | None,
    y_top: float | None,
) -> None:
    """Bin one group's values into its panel and apply the shared y limit, if any."""
    draw_histogram_panel(
        ax,
        data,
        bins=bins,
        log_scale=log_scale,
        xlabel=xlabel,
        ylabel=panel_ylabel,
        title=panel_title,
        bin_edges=bin_edges,
    )
    if y_top is not None and not data.empty:
        ax.set_ylim(0, y_top)


def _draw_grouped_marker(
    ax: Axes,
    *,
    col_val: Any,
    marker_value: float | None,
    marker_label: str,
    marker_on_col_value: str | None,
) -> None:
    """Draw the cutoff marker, only on the panels the caller scoped it to."""
    if marker_value is None:
        return
    if marker_on_col_value is not None and col_val != marker_on_col_value:
        return

    ax.axvline(marker_value, color=FURNITURE_COLOR, linestyle="--", label=marker_label)
    ax.legend(loc="upper right")


@logger.catch(reraise=True)
def render_grouped_histogram_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    value_column: str,
    row_key: str,
    col_key: str | None = None,
    bins: int | np.ndarray = 50,
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

    ``col_key=None`` renders the un-faceted case: one panel per row value, the
    row value as the panel title and the plain ylabel beneath it.
    """
    if col_key is None:
        col_key = _SINGLE_COLUMN
        df = df.assign(**{_SINGLE_COLUMN: ""})
    un_faceted = col_key == _SINGLE_COLUMN

    require_columns(df, [row_key, col_key, value_column], context="grouped histogram input")

    if df.empty:
        logger.warning("No values to plot!")
        return

    grouped = df.groupby([row_key, col_key], sort=True)
    if grouped.ngroups == 0:
        logger.warning("No groups to plot!")
        return

    panel_values = _collect_panel_values(
        grouped, value_column=value_column, log_scale=log_scale, upper_quantile=upper_quantile
    )
    if all(data.empty for data in panel_values.values()):
        logger.warning("All panels are empty after filtering!")
        return

    apply_house_style()

    n_cols = df[col_key].nunique()
    n_rows = df[row_key].nunique()

    # First-appearance order, not sorted: the caller builds the strata in the
    # order they should read (e.g. "All genes" then most-common category first),
    # which alphabetical sorting would scramble.
    row_values = list(dict.fromkeys(df[row_key]))
    col_values = sorted(df[col_key].unique())

    labels = panel_labels(n_rows * n_cols)
    axes = grid_axes(n_rows, n_cols, labels=labels, shape=shape)

    # A row/col combination absent from the data keeps its cell but shows no
    # frame, so the remaining panels stay in their own row and column.
    _hide_absent_cells(axes, panel_values, row_values, col_values, n_cols)

    shared_edges = _shared_log_edges(
        panel_values, share_x_range=share_x_range, log_scale=log_scale, bins=bins
    )
    y_top = _shared_y_top(
        panel_values, share_y_range=share_y_range, shared_edges=shared_edges, bins=bins
    )

    for (row_val, col_val), data in panel_values.items():
        r = row_values.index(row_val)
        c = col_values.index(col_val)
        ax = axes[r * n_cols + c]

        panel_ylabel, panel_title = _grouped_panel_titles(
            row_val, col_val, c, un_faceted=un_faceted, ylabel=ylabel
        )
        _draw_grouped_panel(
            ax,
            data,
            panel_ylabel=panel_ylabel,
            panel_title=panel_title,
            bins=bins,
            log_scale=log_scale,
            xlabel=xlabel,
            bin_edges=shared_edges,
            y_top=y_top,
        )

        logger.info(f"  Panel ({r},{c}): {row_val} {col_val} (n={len(data):,})")
        _draw_grouped_marker(
            ax,
            col_val=col_val,
            marker_value=marker_value,
            marker_label=marker_label,
            marker_on_col_value=marker_on_col_value,
        )

    fit_panels()

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
