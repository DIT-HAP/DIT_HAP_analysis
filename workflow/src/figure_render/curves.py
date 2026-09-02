"""Generic observed-vs-fitted curve grid rendering.

Supersedes ``curve_fitting.py``. The curve model is injected by the caller, so
this module holds no specific model; the pipeline passes
``depletion.curve_model.sigmoid_function``.

Value columns come from the caller rather than a name whitelist. The old
renderer selected them with a literal ['YES0'..'YES4'] test, which silently
picked the wrong subset for any project with a different timepoint count and
raised "x and y must be the same size" from inside matplotlib.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-01
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Callable, Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from loguru import logger

from workflow.src.figures import (
    PanelShape,
    apply_house_style,
    fit_panels,
    grid_axes,
    observed_fitted_colors,
    panel_labels,
    save_dual,
)

from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================
# Updated y-axis range to accommodate LFC values from -10 to 3
# (changed from the original -1.5 to 8.5 to support newer data format)
DEFAULT_YLIM: tuple[float, float] = (-10.5, 3.5)

# Points used to draw the fitted curve smoothly between the observed x values.
FITTED_CURVE_RESOLUTION = 100


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_fitted_curves_figure(
    observed: pd.DataFrame,
    output_stem: Path,
    *,
    x_values: Sequence[float],
    value_columns: Sequence[str],
    model: Callable[..., np.ndarray],
    model_params: Sequence[str],
    annotations: Sequence[str] = (),
    xlabel: str,
    ylabel: str,
    ylim: tuple[float, float] | None = DEFAULT_YLIM,
    n_cols: int = 4,
    shape: PanelShape = PanelShape.WIDE,
    title_column: str | None = None,
) -> None:
    """Render one panel per row: observed points plus the model curve from its parameters.

    Panels default to WIDE: x is a handful of timepoints and y is the fitted
    curve, so the extra width is what separates the observed points.

    Args:
        title_column: Optional column name to use for panel titles instead of the index.
                     Useful for insertion data where titles include gene context.
    """
    logger.info("Rendering fitted curves figure...")

    # Checked before the empty guard so a mismatch is reported even for an empty
    # frame, and before any column access so the message names the real problem.
    if len(x_values) != len(value_columns):
        raise ValueError(
            f"x_values and value_columns must have equal length: "
            f"got {len(x_values)} x values and {len(value_columns)} value columns "
            f"({list(value_columns)})"
        )

    require_columns(
        observed,
        [*value_columns, *model_params, *annotations],
        context="fitted curves input",
    )

    if observed.empty:
        logger.warning("No data to plot!")
        return

    apply_house_style()

    n_panels = len(observed)
    n_rows = (n_panels + n_cols - 1) // n_cols
    logger.info(f"Creating figure with {n_rows}x{n_cols} grid ({n_panels} panels)...")

    observed_color, fitted_color = observed_fitted_colors()

    labels = panel_labels(n_panels)
    axes = grid_axes(n_rows, n_cols, labels=labels, shape=shape)

    x_array = np.asarray(x_values, dtype=float)
    x_smooth = np.linspace(x_array.min(), x_array.max(), FITTED_CURVE_RESOLUTION)

    for ax, (row_index, row) in zip(axes, observed.iterrows(), strict=False):
        y_observed = row[list(value_columns)].to_numpy(dtype=float)
        params = [float(row[name]) for name in model_params]

        ax.scatter(x_array, y_observed, s=30, color=observed_color, alpha=0.8)
        ax.plot(x_smooth, model(x_smooth, *params), color=fitted_color, label="Fitted")

        if annotations:
            text = "\n".join(f"{name}={row[name]:.3f}" for name in annotations)
            ax.text(0.05, 0.95, text, transform=ax.transAxes, verticalalignment="top")

        # Use title_column if provided, otherwise format the index
        if title_column and title_column in row:
            title = str(row[title_column])
        else:
            title = row_index if isinstance(row_index, str) else " ".join(map(str, np.atleast_1d(row_index)))
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        if ylim is not None:
            ax.set_ylim(*ylim)

    for ax in axes[n_panels:]:
        ax.set_visible(False)

    fit_panels()

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
