"""
Generic Plotting
================

General-purpose charts with NO biology-specific concepts (no gene IDs, no
genomics assumptions) — safe to reuse across any project. Ported from
DIT_HAP_pipeline/workflow/src/plot.py (design doc §7: `plotting/generic.py` must
stay domain-agnostic; genomics-aware plots live in gene_level.py instead).

Importing this module applies config/DIT_HAP.mplstyle via plotting.style, so
every caller gets the shared publication look without re-applying it.

Usage
-----
    from workflow.src.plotting.generic import create_scatter_correlation_plot, donut_chart
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Data Processing Imports
import numpy as np
import pandas as pd

# 2. Third-party Imports
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

# 3. Local Imports
# Imported for its side effect: applies config/DIT_HAP.mplstyle once (design doc §7).
from workflow.src.plotting import style as _style  # noqa: F401


# =============================================================================
# SCATTER CORRELATION
# =============================================================================
def create_scatter_correlation_plot(
    x: pd.Series | np.ndarray | list,
    y: pd.Series | np.ndarray | list,
    ax: Axes,
    xscale: None | str = None,
    yscale: None | str = None,
    show_diagonal: bool = True,
    **kwargs,
) -> Axes:
    """Scatter two vectors with a y=x diagonal and a PCC/R²/slope/RMSE stats box."""
    # Drop NaNs, then (for log axes) non-positive values that can't be log-scaled.
    x, y = np.array(x), np.array(y)
    x, y = x[~np.isnan(x) & ~np.isnan(y)], y[~np.isnan(x) & ~np.isnan(y)]

    mask = np.isfinite(x) & np.isfinite(y)
    if xscale == "log":
        mask &= x > 0
    if yscale == "log":
        mask &= y > 0

    x = x[mask]
    y = y[mask]

    ax.scatter(
        x, y,
        alpha=0.5,
        s=10,
        facecolor="none",
        edgecolor="gray",
        rasterized=True,
        **kwargs,
    )

    # Diagonal reference line (y=x), computed from current axis limits.
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    if show_diagonal:
        min_val = min(min(xlim), min(ylim))
        max_val = max(max(xlim), max(ylim))
        ax.plot([min_val, max_val], [min_val, max_val], "k--", alpha=0.8, linewidth=2)

    # Log-transform (for stats) mirrors the visual scale on each axis.
    if xscale == "log":
        ax.set_xscale("log")
        x_for_fitting = np.log10(x)
    else:
        x_for_fitting = x

    if yscale == "log":
        ax.set_yscale("log")
        y_for_fitting = np.log10(y)
    else:
        y_for_fitting = y

    # Pearson correlation + R².
    pcc = np.corrcoef(x_for_fitting, y_for_fitting)[0, 1]
    r_squared = pcc**2
    try:
        slope, intercept = np.polyfit(x_for_fitting, y_for_fitting, 1)
        y_pred = intercept + slope * x_for_fitting
        rmse = np.sqrt(np.mean((y_for_fitting - y_pred) ** 2))
    except Exception as e:  # noqa: BLE001 — degenerate input (e.g. <2 points) shouldn't abort the figure
        print("Error in linear regression:", e)
        slope, intercept, rmse = np.nan, np.nan, np.nan

    stats_text = [
        f"Data points: {len(x):,}",
        f"PCC: {pcc:.4f}",
        f"R²: {r_squared:.4f}",
        f"Slope: {slope:.4f}",
        f"Intercept: {intercept:.4f}",
        f"RMSE: {rmse:.4f}",
    ]
    ax.text(0.02, 0.98, "\n".join(stats_text), transform=ax.transAxes, verticalalignment="top")

    return ax


# =============================================================================
# DONUT CHART
# =============================================================================
def donut_chart(
    values: list[int],
    labels: list[str],
    colors: list[str],
    center_text: str = "",
    ax: Axes | None = None,
) -> Axes | Figure:
    """Donut (ring) chart with per-wedge percent+count labels and centered text."""
    return_ax = True
    if ax is None:
        fig, ax = plt.subplots()
        return_ax = False

    ax.pie(
        values,
        colors=colors,
        autopct=lambda pct: f"{pct:.1f}%\n({int(round(pct / 100 * sum(values))):,})",
        startangle=90,
        pctdistance=0.75,
        wedgeprops=dict(width=0.5, edgecolor="white"),
        textprops={"fontsize": 22, "weight": "bold"},
    )
    ax.text(0, 0, center_text, ha="center", va="center", fontsize=26, fontweight="bold")
    ax.axis("equal")

    return ax if return_ax else fig
