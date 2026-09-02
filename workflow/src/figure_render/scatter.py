"""Generic scatter and regression figure rendering.

Supersedes ``correlation.py`` (PBL vs PBR) and ``orientation.py`` (plus vs minus
strand), which rendered the same grouped log-log regression grid for different
columns. Column names and axis labels are supplied by the caller.

``ScatterPanel.scale`` owns the log handling: ``"log"`` takes raw columns, puts
both axes on true log scale with exponential major labels plus minor ticks, and
computes r, P and density in log10 space over strictly positive pairs. Callers
therefore pass raw values, not pre-logged columns. ``"symlog"`` is the
display-only variant for data containing zeros; its statistics stay in linear
space.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-01
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from loguru import logger
from matplotlib.axes import Axes
from matplotlib.collections import PathCollection
from matplotlib.colors import Normalize
from scipy.stats import pearsonr

from figures import (
    FURNITURE_COLOR,
    apply_house_style,
    apply_log_scale,
    apply_symlog_scale,
    fit_panels,
    grid_axes,
    panel_labels,
    save_dual,
)

from ._point_density import point_density
from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================
# Verified to render ~50k points per panel legibly: density stays visible instead
# of collapsing into a solid block. cns.scatterplot always forces edgecolor=None
# internally, so a filled gray marker (not regplot's hollow-circle look) is the
# closest equivalent that doesn't collide with that forced kwarg.
REGRESSION_PANEL_SCATTER_KWS: dict[str, object] = {
    "s": 3,
    "color": FURNITURE_COLOR,
    "alpha": 0.15,
    "rasterized": True,
}

# Verified for ~5k points with categorical colour-coding. Deliberately NOT shared
# with REGRESSION_PANEL_SCATTER_KWS: cns.scatterplot always passes edgecolor=None
# to seaborn internally, so supplying 'edgecolor' here raises on the duplicate.
SCATTERPLOT_KWS: dict[str, object] = {
    "s": 6,
    "alpha": 0.4,
    "rasterized": True,
}

# Density panels bypass cns.scatterplot for ax.scatter (cns.scatterplot exposes
# only categorical hue and forces edgecolor=None, so a per-point colour array
# cannot be threaded through it). Higher alpha than
# REGRESSION_PANEL_SCATTER_KWS: at 0.15 the colour signal washes out, and
# overplotting is what the colour already encodes. No 'color'/'c' key -- those
# collide with the density array.
DENSITY_SCATTER_KWS: dict[str, Any] = {
    "s": 1.5,
    "alpha": 0.6,
    "linewidths": 0,
    "rasterized": True,
}

# Axes-fraction (x0, y0, width, height) for the inset colorbar: bottom-right,
# clear of the top-left n/r/P annotation. An inset axes contributes nothing to
# the enclosing grid cell, so it cannot reflow the grid the way
# fig.colorbar(ax=...) would by taking space from the panel itself.
DENSITY_CBAR_BOUNDS = (0.82, 0.02, 0.035, 0.30)
DENSITY_CBAR_LABEL = "Density"

# Panels per row in render_scatter_grid_figure. See the comment at the n_cols
# assignment for the measurement behind 2.
MAX_GRID_COLUMNS = 2

# Right-hand strip reserved for the shared legend, as a tight_layout rect and the
# matching axes-fraction anchor. 0.86 leaves the legend clear of the panels while
# costing the grid little width.
_LEGEND_STRIP_LEFT = 0.86
LEGEND_LAYOUT_RECT = (0, 0, _LEGEND_STRIP_LEFT, 1)
LEGEND_ANCHOR = (_LEGEND_STRIP_LEFT, 0.5)

# Layout points between the bar's top edge and its label.
DENSITY_CBAR_LABEL_PAD = 2

# Colorbar-axes fraction the "Density" title starts at. Negative: the word is
# several bar-widths wide, so it has to begin left of the bar to sit over it.
DENSITY_CBAR_TITLE_X = -1.5

# Legend handles inherit the scatter's small marker and low alpha, so they are
# scaled up and forced opaque or the swatches misrepresent the plotted colours.
LEGEND_MARKER_SCALE = 2


@dataclass(kw_only=True, slots=True, frozen=True)
class ScatterPanel:
    """One panel of a scatter grid: which columns, how to label, which reference line."""

    x: str
    y: str
    xlabel: str
    ylabel: str
    title: str
    reference: Literal["identity", "zero", "unit_identity", "none"] = "none"
    scale: Literal["linear", "log", "symlog"] = "linear"
    show_stats: bool = False
    density: bool = False


# =============================================================================
# CORE LOGIC
# =============================================================================
def _has_hue(df: pd.DataFrame, hue: str | None) -> bool:
    """Report whether hue colouring applies, matching render_scatter_panel's own test."""
    return hue is not None and hue in df.columns



def _annotate_fit_stats(ax: Axes, df: pd.DataFrame, *, x: str, y: str) -> None:
    """Annotate n, r, and P at regplot's own text position/style, so toggling show_stats doesn't reflow the panel."""
    valid = df[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(valid) < 2:
        logger.warning(f"Not enough finite pairs ({len(valid)}) to compute r/P for {x!r} vs {y!r}")
        return

    r, p = pearsonr(valid[x], valid[y])

    if p == 0.0 or p < 1e-9:
        p_text = rf"$P < 1\mathrm{{E}}{{-}}9$"
    else:
        # Format as mantissa E exponent
        exponent = int(np.floor(np.log10(p)))
        mantissa = p / (10 ** exponent)
        p_text = rf"$P={mantissa:.2f}\mathrm{{E}}{{{exponent:+03d}}}$"

    ax.text(
        0.05,
        0.95,
        rf"$n={len(valid)}$" "\n" rf"$r={r:.2f}$" "\n" + p_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
    )


def _draw_reference_line(ax: Axes, df: pd.DataFrame, panel: ScatterPanel) -> None:
    """Draw the panel's reference guide line, dashed gray to read as a guide rather than a fit.

    Fitted/regression lines (drawn elsewhere via cns.regplot or curves.py) stay
    solid so the two kinds of line remain visually distinct.
    """
    positive = df[(df[panel.x] > 0) & (df[panel.y] > 0)] if panel.scale == "log" else df

    match panel.reference:
        case "identity":
            min_val = min(positive[panel.x].min(), positive[panel.y].min())
            max_val = max(positive[panel.x].max(), positive[panel.y].max())
            ax.plot([min_val, max_val], [min_val, max_val], color=FURNITURE_COLOR, linestyle="--")
        case "unit_identity":
            ax.plot([0, 1], [0, 1], color=FURNITURE_COLOR, linestyle="--")
        case "zero":
            ax.axhline(0, color=FURNITURE_COLOR, linestyle="--")
        case "none":
            pass


def _apply_shared_square_limits(
    axes: Sequence[Axes], df: pd.DataFrame, *, x: str, y: str, scale: Literal["linear", "log", "symlog"] = "linear",
) -> None:
    """Put every axes on one common x=y range spanning all finite data.

    Alignment of the axes boxes is geometric (see ``figures.grid_axes``); this
    makes the *contents* comparable too, which is what lets a reader read one
    identity line across the whole grid. A common range on both axes also keeps
    the identity reference at 45 degrees in every panel.
    """
    finite = df[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
    if scale == "log":
        finite = finite[(finite[x] > 0) & (finite[y] > 0)]
    if finite.empty:
        logger.warning("No finite positive values; leaving panel limits autoscaled")
        return

    low = float(min(finite[x].min(), finite[y].min()))
    high = float(max(finite[x].max(), finite[y].max()))
    if not np.isfinite(low) or not np.isfinite(high) or low == high:
        logger.warning(f"Degenerate data range [{low}, {high}]; leaving panel limits autoscaled")
        return

    if scale == "log":
        margin_factor = (high / low) ** 0.04
        for ax in axes:
            ax.set_xlim(low / margin_factor, high * margin_factor)
            ax.set_ylim(low / margin_factor, high * margin_factor)
    else:
        margin = (high - low) * 0.04
        for ax in axes:
            ax.set_xlim(low - margin, high + margin)
            ax.set_ylim(low - margin, high + margin)


def _add_density_colorbar(ax: Axes, mappable: PathCollection) -> None:
    """Inset a thin vertical colorbar labelled low/high inside the panel's bottom-right corner."""
    cax = ax.inset_axes(DENSITY_CBAR_BOUNDS)
    colorbar = ax.figure.colorbar(mappable, cax=cax, ticks=[0, 1])

    # Only the ordering is meaningful: density is renormalised per panel, so
    # numeric ticks would invite comparison between panels with different n.
    labels = colorbar.ax.set_yticklabels(["low", "high"])

    # End ticks anchor their label on the tick centre, pushing half the text
    # past the bar. Aligning each label inwards keeps both inside the panel.
    labels[0].set_verticalalignment("bottom")
    labels[-1].set_verticalalignment("top")

    # length=0 drops the tick marks but keeps their labels.
    colorbar.ax.tick_params(length=0, pad=1)

    # Sits above the bar, sized and right-aligned to match the low/high ticks: the
    # bar is far narrower than the word, so a centred title would spill past the
    # panel's right edge, and an axis label would read as a full-size axis title.
    title = colorbar.ax.set_title(
        DENSITY_CBAR_LABEL,
        loc="center",
        # fontsize=plt.rcParams["ytick.labelsize"],
        # pad=DENSITY_CBAR_LABEL_PAD,
    )
    # title.set_position((DENSITY_CBAR_TITLE_X, title.get_position()[1]))
    colorbar.outline.set_linewidth(cns.settings.axes_linewidth)


def _draw_density_scatter(ax: Axes, df: pd.DataFrame, panel: ScatterPanel) -> bool:
    """Colour the panel's points by 2D density, returning False when the cloud is too degenerate.

    Deliberately ignores the caller's ``scatter_kws``: those are tuned for the
    cns.scatterplot path and carry a flat ``color`` that collides with the
    per-point colour array.
    """
    x_vals = np.log10(df[panel.x].to_numpy()) if panel.scale == "log" else df[panel.x].to_numpy()
    y_vals = np.log10(df[panel.y].to_numpy()) if panel.scale == "log" else df[panel.y].to_numpy()
    density = point_density(x_vals, y_vals)

    if density is None:
        return False

    finite = np.isfinite(density)
    if panel.scale == "log":
        finite = finite & (df[panel.x] > 0).to_numpy() & (df[panel.y] > 0).to_numpy()

    # Ascending so the crowded points draw last and are not buried under the
    # sparse ones.
    order = np.argsort(density[finite])
    x = df[panel.x].to_numpy()[finite][order]
    y = df[panel.y].to_numpy()[finite][order]

    # Rescaled to [0, 1] so the colorbar's low/high ticks sit at the bar ends.
    # A bare Normalize() autoscales to the raw density range instead, leaving
    # ticks at 0 and 1 out of range and clipped onto the edges.
    ranked = density[finite][order]
    span = np.ptp(ranked)
    scaled = (ranked - ranked.min()) / span if span > 0 else np.zeros_like(ranked)

    # Colormap comes from rcParams['image.cmap'], which apply_house_style() points
    # at the house sequential map.
    mappable = ax.scatter(
        x, y, c=scaled, norm=Normalize(vmin=0, vmax=1), **DENSITY_SCATTER_KWS
    )
    _add_density_colorbar(ax, mappable)

    return True


def render_scatter_panel(
    ax: Axes,
    df: pd.DataFrame,
    panel: ScatterPanel,
    *,
    hue: str | None = None,
    hue_order: Sequence[str] | None = None,
    scatter_kws: Mapping[str, object] | None = None,
    show_legend: bool = True,
) -> None:
    """Draw one ScatterPanel (or an empty-data placeholder) onto ax.

    The single-axes primitive behind both ``render_scatter_grid_figure`` and
    ``render_grouped_regression_figure``. Any other grid orchestrator that
    wants one cns.scatterplot panel with this project's hue/reference-line/
    scale/stats handling can call this directly instead of duplicating
    that sequence.

    ``panel.scale="log"`` applies true log axes with minor ticks and computes
    stats/density in log10 space on strictly positive values; caller passes
    **raw** columns (e.g. ``pbl``), not pre-logged columns.

    ``panel.show_stats=True`` annotates n, r, and P (computed over all of df,
    ignoring hue) in the same text position/style cnsplots' own regplot used.

    ``panel.density=True`` colours points by 2D kernel density with an inset
    colorbar. It is ignored when ``hue`` applies: a categorical hue and a
    continuous density both claim the marker colour.

    ``show_legend=False`` suppresses the per-panel legend; callers creating
    multi-panel grids should pass False and create a shared figure-level legend.
    """
    ax.set_xlabel(panel.xlabel)
    ax.set_ylabel(panel.ylabel)
    ax.set_title(panel.title)

    if df.empty:
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
        return

    plot_df = df[(df[panel.x] > 0) & (df[panel.y] > 0)].copy() if panel.scale == "log" else df.copy()
    if plot_df.empty:
        ax.text(0.5, 0.5, "No positive data", ha="center", va="center", transform=ax.transAxes)
        return

    kws = dict(SCATTERPLOT_KWS if scatter_kws is None else scatter_kws)

    if hue is not None and hue in plot_df.columns:
        if panel.density:
            logger.warning(
                f"Panel '{panel.title}': hue {hue!r} and density both set; "
                f"colouring by hue and skipping density"
            )
        # Levels absent from the data must be dropped, or a project lacking
        # one of them raises inside seaborn.
        observed = set(plot_df[hue].unique())
        order = [level for level in (hue_order or sorted(observed)) if level in observed]
        cns.scatterplot(data=plot_df, x=panel.x, y=panel.y, hue=hue, hue_order=order, ax=ax, **kws)
        # seaborn draws its own legend, so suppressing ours is not enough. Removing
        # the drawn artist keeps get_legend_handles_labels() working for a shared
        # legend, which passing legend=False to seaborn would leave empty.
        if show_legend:
            ax.legend(loc="best")
        elif ax.get_legend() is not None:
            ax.get_legend().remove()
    else:
        density_drawn = panel.density and _draw_density_scatter(ax, plot_df, panel)
        if panel.density and not density_drawn:
            logger.warning(f"Panel '{panel.title}': falling back to flat colour, no density estimate")
        if not density_drawn:
            cns.scatterplot(data=plot_df, x=panel.x, y=panel.y, ax=ax, **kws)

    ax.set_xlabel(panel.xlabel)
    ax.set_ylabel(panel.ylabel)
    ax.set_title(panel.title)

    _draw_reference_line(ax, plot_df, panel)

    if panel.show_stats:
        if panel.scale == "log":
            log_df = plot_df.copy()
            log_df["_log10_x"] = np.log10(plot_df[panel.x])
            log_df["_log10_y"] = np.log10(plot_df[panel.y])
            _annotate_fit_stats(ax, log_df, x="_log10_x", y="_log10_y")
        else:
            _annotate_fit_stats(ax, plot_df, x=panel.x, y=panel.y)

    if panel.scale == "log":
        apply_log_scale(ax)
    elif panel.scale == "symlog":
        apply_symlog_scale(ax)


@logger.catch(reraise=True)
def render_grouped_regression_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    x: str,
    y: str,
    xlabel: str,
    ylabel: str,
    row_key: str = "sample",
    col_key: str = "timepoint",
    show_stats: bool = True,
    density: bool = False,
    share_limits: bool = True,
    scale: Literal["linear", "log", "symlog"] = "linear",
    scatter_kws: Mapping[str, object] | None = None,
) -> None:
    """Render a row_key x col_key grid of square regression panels for the x/y columns.

    show_stats, density and scale are set on every panel in the grid: this
    renderer builds its own ScatterPanel specs, so they cannot be varied per
    panel here.

    ``scale="log"`` takes raw (unlogged) x/y columns and puts both axes on true
    log scale, computing r/P and density in log10 space over strictly positive
    pairs.

    ``share_limits=True`` puts every panel on one common x=y range, so panels
    are comparable by eye as well as aligned; pass False to let each panel
    autoscale.
    """
    logger.info("Rendering grouped regression figure...")

    # Validated even when df is empty, so a typo'd column name is caught rather
    # than silently returning with no artifact and no error.
    require_columns(df, [row_key, col_key, x, y], context="scatter input")

    if df.empty:
        logger.warning("No data to plot!")
        return

    apply_house_style()

    grouped = df.groupby([row_key, col_key], sort=True)
    n_panels = grouped.ngroups

    if n_panels == 0:
        logger.warning("No groups to plot!")
        return

    logger.info(f"Creating figure with {n_panels} panels...")

    # One row per row_key value and one column per col_key value, so a row stays
    # on one line instead of wrapping mid-group.
    row_values = sorted(df[row_key].unique())
    col_values = sorted(df[col_key].unique())
    n_rows, n_cols = len(row_values), len(col_values)

    labels = panel_labels(n_rows * n_cols)
    axes = grid_axes(n_rows, n_cols, labels=labels)

    groups = dict(iter(grouped))

    for row_index, row_value in enumerate(row_values):
        for col_index, col_value in enumerate(col_values):
            panel_index = row_index * n_cols + col_index
            label = labels[panel_index]
            ax = axes[panel_index]

            # An absent row/col combination still gets its cell, so the grid keeps
            # one column per col_key value instead of shifting later panels left.
            group_df = groups.get((row_value, col_value), df.iloc[:0])

            # The row_key value rides in the first column's ylabel rather than the
            # title: a full "{row} {col}" title is wider than the axes and would
            # overrun the neighbouring panel on this fixed pitch.
            panel_ylabel = f"{row_value}\n{ylabel}" if col_index == 0 else ylabel

            if group_df.empty:
                logger.warning(f"  Panel {label}: {row_value} {col_value} has no valid data")
            else:
                logger.info(f"  Panel {label}: {row_value} {col_value} (n={len(group_df)})")

            panel = ScatterPanel(
                x=x, y=y, xlabel=xlabel, ylabel=panel_ylabel, title=str(col_value),
                reference="identity", show_stats=show_stats, density=density, scale=scale,
            )
            render_scatter_panel(
                ax, group_df, panel,
                scatter_kws=REGRESSION_PANEL_SCATTER_KWS if scatter_kws is None else scatter_kws,
            )

    if share_limits:
        _apply_shared_square_limits(axes, df, x=x, y=y, scale=scale)

    fit_panels()

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")


@logger.catch(reraise=True)
def render_scatter_grid_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    panels: Sequence[ScatterPanel],
    hue: str | None = None,
    hue_order: Sequence[str] | None = None,
    scatter_kws: Mapping[str, object] | None = None,
) -> None:
    """Render one scatter panel per ScatterPanel spec, with optional hue colouring."""
    logger.info(f"Rendering scatter grid with {len(panels)} panels...")

    if not panels:
        logger.warning("No panels requested!")
        return

    # Validated even when df is empty, so a typo'd column name is caught rather
    # than silently rendering a placeholder panel.
    for panel in panels:
        require_columns(df, [panel.x, panel.y], context=f"scatter panel '{panel.title}'")

    apply_house_style()

    # Same plot type in every cell, so the panels belong in an aligned grid.
    # Building one from repeated multipanel.panel() calls would leave each cell
    # offset by the rendered width of its own y tick labels.
    #
    # 2-across, not 4: at 4 columns a 540 px page leaves each axes 182 device px
    # while these titles render 204-286 px, so every title overruns its own panel
    # and collides with the next panel's letter. 2 columns give 378 px, which
    # clears the widest measured title.
    n_cols = min(len(panels), MAX_GRID_COLUMNS)
    n_rows = (len(panels) + n_cols - 1) // n_cols

    labels = panel_labels(len(panels))
    axes = grid_axes(n_rows, n_cols, labels=labels)

    for label, panel, ax in zip(labels, panels, axes, strict=True):
        if df.empty:
            logger.warning(f"  Panel {label}: {panel.title} has no valid data")
        else:
            logger.info(f"  Panel {label}: {panel.title} (n={len(df)})")

        render_scatter_panel(ax, df, panel, hue=hue, hue_order=hue_order, scatter_kws=scatter_kws, show_legend=False)

    for ax in axes[len(panels):]:
        ax.set_visible(False)

    fig = plt.gcf()

    # fit_panels before the legend: it measures each panel's box to ensure they
    # match their PanelShape. Adding the legend first would include it in those
    # measurements.
    fit_panels(rect=LEGEND_LAYOUT_RECT if _has_hue(df, hue) else None)

    # One figure-level legend outside the axes. Per-panel legends repeat the same
    # hue levels in every cell and, at loc="best", land on the data or the stats
    # annotation, since "best" only searches within one axes.
    if _has_hue(df, hue):
        handles, legend_labels = axes[0].get_legend_handles_labels()
        if handles:
            legend = fig.legend(
                handles, legend_labels,
                loc="center left", bbox_to_anchor=LEGEND_ANCHOR,
                frameon=False, title=hue, markerscale=LEGEND_MARKER_SCALE,
                alignment = "left"
            )
            # Left-align and position the title over the entry labels rather than
            # centred over the whole box, so it reads as their heading.
            # title = legend.get_title()
            # title.set_horizontalalignment("left")
            # title.set_position((0, 0))  # Align with legend entries
            for handle in legend.legend_handles:
                handle.set_alpha(1.0)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
