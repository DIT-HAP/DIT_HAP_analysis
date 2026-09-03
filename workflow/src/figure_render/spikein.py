"""
Spike-In Panel Renderer (cnsplots-adapted)
===========================================

cnsplots-style renderer for spike-in dilution linearity panel. Replaces the
traditional matplotlib implementation in workflow/src/pcr_qc/core.py with a
figure_render-compatible approach.

This is a domain-specific renderer for spike-in QC that cannot use a generic
figure_render plotter (no built-in spike-in plot type), so it follows cnsplots
conventions directly.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-02
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from loguru import logger

# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def render_spikein_panel(
    ax: plt.Axes,
    spikein: pd.DataFrame,
    *,
    marker_size: float = 150,
    marker_linewidth: float = 1.5,
    marker_alpha: float = 0.9,
    fit_line_color: str = "black",
    fit_line_style: str = "--",
    fit_line_width: float = 2.5,
    fit_line_alpha: float = 0.7,
    show_stats: bool = True,
    stats_position: tuple[float, float] = (0.05, 0.95),
    legend_location: str = "lower right",
) -> None:
    """Render spike-in dilution linearity panel with cnsplots style.

    Draws a scatter plot of log2(relative dilution ratio) vs log2(relative read
    ratio) for each spike-in insertion site, plus a linear regression fit line.
    Excludes Spikein0 (zero-dilution reference) from the fit.

    Parameters
    ----------
    ax : plt.Axes
        Target axes (already styled by cnsplots)
    spikein : pd.DataFrame
        Spike-in stats with columns: Sample, Name, Relative_Dilution_Ratio,
        Relative_Read_Ratio
    marker_size : float, default 150
        Scatter marker size in points²
    marker_linewidth : float, default 1.5
        Scatter marker edge width
    marker_alpha : float, default 0.9
        Scatter marker alpha
    fit_line_color : str, default "black"
        Linear fit line color
    fit_line_style : str, default "--"
        Linear fit line style
    fit_line_width : float, default 2.5
        Linear fit line width
    fit_line_alpha : float, default 0.7
        Linear fit line alpha
    show_stats : bool, default True
        Whether to show fit statistics (PCC, R², Slope, Intercept)
    stats_position : tuple[float, float], default (0.05, 0.95)
        Statistics text position in axes coordinates (x, y)
    legend_location : str, default "lower right"
        Legend location

    Notes
    -----
    - This function assumes ax is already styled by cnsplots (via apply_house_style)
    - Uses categorical colors from the current color cycle
    - Marker colors cycle through available palette colors
    - Does NOT call apply_house_style itself (caller's responsibility)
    """
    # Exclude Spikein0 (zero-dilution reference) from fit, matching original logic
    spikein_filtered = spikein.query("Sample != 'Spikein0'")

    # Get current color cycle from rcParams (set by cnsplots)
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    # Scatter plot for each spike-in insertion site
    for idx, (name, sub) in enumerate(spikein_filtered.groupby("Name")):
        color = color_cycle[idx % len(color_cycle)]
        ax.scatter(
            sub["Relative_Dilution_Ratio"],
            sub["Relative_Read_Ratio"],
            label=name,
            facecolor="none",
            edgecolor=color,
            s=marker_size,
            linewidth=marker_linewidth,
            alpha=marker_alpha,
        )

    # Linear regression fit
    slope, intercept, r_value, p_value, std_err = linregress(
        spikein_filtered["Relative_Dilution_Ratio"],
        spikein_filtered["Relative_Read_Ratio"],
    )

    # Draw fit line
    x_range = spikein_filtered["Relative_Dilution_Ratio"]
    x_min, x_max = x_range.min(), x_range.max()
    line_x = np.array([x_min - 1, x_max + 1])  # Extend slightly beyond data
    line_y = slope * line_x + intercept

    ax.plot(
        line_x,
        line_y,
        color=fit_line_color,
        linestyle=fit_line_style,
        linewidth=fit_line_width,
        alpha=fit_line_alpha,
        zorder=1,  # Behind scatter points
    )

    # Set axis labels (mathematical notation)
    ax.set_xlabel(r"$\log_2$(relative dilution ratio)")
    ax.set_ylabel(r"$\log_2$(relative read ratio)")

    # Show fit statistics
    if show_stats:
        stats_text = (
            f"PCC = {r_value:.2f}\n"
            f"R² = {r_value**2:.2f}\n"
            f"Slope = {slope:.2f}\n"
            f"Intercept = {intercept:.2f}"
        )
        ax.text(
            stats_position[0],
            stats_position[1],
            stats_text,
            transform=ax.transAxes,
            ha="left",
            va="top",
            # Do not specify fontsize - inherit from cnsplots rcParams
        )

    # Legend (frameon controlled by cnsplots rcParams)
    ax.legend(loc=legend_location)

    # Set reasonable tick positions if data range is known
    # (cnsplots handles tick formatting)
    if x_range.min() >= -10 and x_range.max() <= 0:
        ax.set_xticks(np.arange(-10, 1, 2))

    logger.debug(f"Rendered spike-in panel: {len(spikein_filtered)} points, R²={r_value**2:.3f}")
