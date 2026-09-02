"""Generic MA plot rendering.

Supersedes ``ma_plot.py`` (wide baseMean/LFC tables, one column per timepoint)
and ``ma_plot_replicates.py`` (one long table with padj-driven point colours).
Callers reshape either form into a sequence of (title, abundance, effect)
triples, so this module needs no knowledge of the input layout.

The two pipeline branches are selected by ``has_replicates`` and never run
together; both must keep working.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-01
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import sys
from collections.abc import Sequence
from enum import StrEnum
from pathlib import Path

# Add workflow/src to path for imports
_SCRIPT_DIR = Path(__file__).resolve().parent
_SRC_DIR = _SCRIPT_DIR.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger

from figures import (  # noqa: E402
    FURNITURE_COLOR,
    PanelShape,
    apply_house_style,
    apply_log_scale,
    fit_panels,
    grid_axes,
    panel_labels,
    save_dual,
)


# =============================================================================
# CONSTANTS
# =============================================================================
class Orientation(StrEnum):
    """Panel arrangement and axis assignment, which move together.

    VERTICAL stacks n panels in one column with the effect on x; HORIZONTAL lays
    them in one row with the abundance on x.
    """

    VERTICAL = "vertical"
    HORIZONTAL = "horizontal"


NONSIGNIFICANT_COLOR = FURNITURE_COLOR

# Rasterize: tens of thousands of points per panel would bloat the vector PDF.
MA_SCATTER_KWS: dict[str, object] = {
    "s": 1.5,
    "alpha": 0.4,
    "linewidths": 0,
    "rasterized": True,
}


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_ma_figure(
    panels: Sequence[tuple[str, pd.Series, pd.Series]],
    output_stem: Path,
    *,
    abundance_label: str,
    effect_label: str,
    title_prefix: str,
    orientation: Orientation = Orientation.VERTICAL,
    point_colors: Sequence[pd.Series | str] | None = None,
    share_axes: bool = True,
    null_effect: float = 0.0,
) -> None:
    """Render one square MA panel per (title, abundance, effect) triple.

    All panels sit in a single line -- one column under VERTICAL, one row under
    HORIZONTAL -- so the panel count alone sets the grid.
    """
    logger.info(f"Rendering MA plot figure ({orientation})...")

    if not panels:
        logger.warning("No data to plot!")
        return

    if point_colors is not None and len(point_colors) != len(panels):
        raise ValueError(
            f"point_colors must have one entry per panel: "
            f"got {len(point_colors)} colours for {len(panels)} panels"
        )

    apply_house_style()

    n_panels = len(panels)
    logger.info(f"Creating figure with {n_panels} panels...")

    stacked = orientation is Orientation.VERTICAL
    n_rows, n_cols = (n_panels, 1) if stacked else (1, n_panels)

    labels = panel_labels(n_panels)
    axes = grid_axes(
        n_rows, n_cols, labels=labels, shape=PanelShape.SQUARE,
        share_x=share_axes, share_y=share_axes,
    )

    colors = point_colors if point_colors is not None else [NONSIGNIFICANT_COLOR] * n_panels

    for ax, label, (title, abundance, effect), color in zip(
        axes, labels, panels, colors, strict=True
    ):
        logger.info(f"  Panel {label}: {title} (n={len(effect)})")

        match orientation:
            case Orientation.VERTICAL:
                # Effect on x, abundance (log) on y: reference line is vertical.
                ax.scatter(effect, abundance, c=color, **MA_SCATTER_KWS)
                apply_log_scale(ax, x=False, y=True)
                ax.axvline(null_effect, color=FURNITURE_COLOR, linestyle="--", zorder=3)
                ax.set_xlabel(effect_label)
                ax.set_ylabel(abundance_label)
            case Orientation.HORIZONTAL:
                # Abundance (log) on x, effect on y: reference line is horizontal.
                ax.scatter(abundance, effect, c=color, **MA_SCATTER_KWS)
                apply_log_scale(ax, x=True, y=False)
                ax.axhline(null_effect, color=FURNITURE_COLOR, linestyle="--", zorder=3)
                ax.set_xlabel(abundance_label)
                ax.set_ylabel(effect_label)

        ax.set_title(f"{title_prefix} - {title}")

    fit_panels()

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
