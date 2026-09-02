"""Generic composition figure rendering: a bar overview plus per-category donuts.

Supersedes ``coverage.py``, which showed gene coverage percentage per viability
category alongside covered/not-covered donuts. Column names and labels are
supplied by the caller.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger

from workflow.src.figures import (
    JOURNAL_HEIGHT_PX,
    JOURNAL_WIDTH_PX,
    PANEL_SIZES_PX,
    PanelShape,
    apply_house_style,
    panel_labels,
    save_dual,
)

from ._schema import require_columns
from .donut import draw_donut_panel

# =============================================================================
# CONSTANTS
# =============================================================================
# This figure is genuinely heterogeneous — one bar panel plus a block of
# repeated donuts — so it uses multipanel for the two blocks and a GridSpec
# inside the donut host. The donuts must align with each other, which repeated
# multipanel.panel() calls cannot deliver: each panel's origin is offset by the
# rendered width of its own y decorations.
#
# Bar and donut cells take explicit px dimensions. They could be expressed as
# PanelShape values, but the bar's pad_left is tuned for the category names
# (longer than the usual ylabel), so it stays as explicit px for now.
BAR_PANEL_WIDTH_PX, BAR_PANEL_HEIGHT_PX = PANEL_SIZES_PX[PanelShape.SQUARE]
BAR_PANEL_PAD_LEFT_PX = 45  # gap for the category names, as in the gallery
BAR_PANEL_MARGIN_RIGHT_PX = 8  # match the donut host's right margin
MULTIPANEL_MAX_WIDTH_PX = 250
DONUT_PANEL_WIDTH_PX, DONUT_PANEL_HEIGHT_PX = PANEL_SIZES_PX[PanelShape.SQUARE]
DONUTS_PER_ROW = 2
DONUT_PANEL_MARGIN_RIGHT_PX = 8

# multipanel measures each panel's *top* decorations (title, panel label) and
# reserves layout space for them, but nothing below the axes: x tick labels,
# xlabel and a legend="bottom" all render outside the panel's total height. With
# the 10 px default the bar's xlabel landed on top of the donut panel labels.
# These reserve the measured overhang (~20 layout px for the bar's ticks +
# xlabel) plus ROW_GAP_PX of air, so the blocks read as separate.
ROW_GAP_PX = 14
BAR_PANEL_MARGIN_BOTTOM_PX = 20 + ROW_GAP_PX
DONUT_PANEL_MARGIN_BOTTOM_PX = 22 + ROW_GAP_PX

# Fractions of one grid cell left between donut cells. A donut's legend sits
# under its ring and its panel label to the upper left, and neither is inside
# the axes box, so the gaps carry that overhang.
DONUT_GRID_WSPACE = 0.35
DONUT_GRID_HSPACE = 0.52

# The host axes is layout scaffolding, not a panel: it is invisible and holds no
# data. It cannot be removed after subdividing — multipanel's draw handler still
# reaches for it and fails — so it is tagged instead, letting callers and tests
# tell scaffolding from panels.
DONUT_GRID_HOST_GID = "donut-grid-host"

# Donut-hole text: at the house legend size the hole fits ~15 chars per line.
# Category names longer than that are hyphen-wrapped before the n-count line.
CENTER_LINE_MAX_CHARS = 15


def _wrap_center_text(center_text: str) -> str:
    """Hyphen-wrap ``center_text`` so no line exceeds CENTER_LINE_MAX_CHARS."""
    if len(center_text) <= CENTER_LINE_MAX_CHARS:
        return center_text
    lines: list[str] = []
    current = ""
    for part in center_text.split("-"):
        candidate = f"{current}-{part}" if current else part
        if len(candidate) > CENTER_LINE_MAX_CHARS and current:
            lines.append(current)
            current = part
        else:
            current = candidate
    if current:
        lines.append(current)
    return "\n".join(lines)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_composition_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    category_column: str,
    percentage_column: str,
    part_column: str,
    whole_column: str,
    part_label: str,
    whole_label: str,
    xlabel: str,
    ylabel: str,
    title: str,
) -> None:
    """Render a horizontal percentage bar chart plus one part/whole donut each.

    Categories run down the y axis so long names stay readable; percentage
    labels sit at the bar ends. Each donut carries its category name in the
    hole instead of a title above, saving the vertical title band.
    """
    logger.info("Rendering composition figure...")

    require_columns(
        df,
        [category_column, percentage_column, part_column, whole_column],
        context="composition input",
    )

    apply_house_style()

    if df.empty:
        logger.warning("No composition data to plot!")
        cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
        multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)
        ax = multipanel.panel(label="A")
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
        save_dual(output_stem)
        return

    n_categories = len(df)
    logger.info(f"Creating figure with {1 + n_categories} panels ({n_categories} categories)...")

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=MULTIPANEL_MAX_WIDTH_PX)
    labels = panel_labels(1 + n_categories)

    n_donut_rows = (n_categories + DONUTS_PER_ROW - 1) // DONUTS_PER_ROW

    # Panel A: horizontal percentage bars, categories down the y axis so long
    # names stay readable; percentages annotate the bar ends.
    ax_bar = multipanel.panel(
        label=labels[0],
        width=BAR_PANEL_WIDTH_PX,
        height=BAR_PANEL_HEIGHT_PX,
        pad_left=BAR_PANEL_PAD_LEFT_PX,
        margin_right=BAR_PANEL_MARGIN_RIGHT_PX,
        margin_bottom=BAR_PANEL_MARGIN_BOTTOM_PX,
    )

    # An empty host rectangle whose coordinates become the donut grid's bounds.
    # Sized for the whole block so the grid inherits one aligned pitch.
    # No label here — the first donut cell will carry labels[1].
    host = multipanel.panel(
        label="",
        width=DONUT_PANEL_WIDTH_PX * DONUTS_PER_ROW,
        height=DONUT_PANEL_HEIGHT_PX * n_donut_rows,
        margin_right=DONUT_PANEL_MARGIN_RIGHT_PX,
        margin_bottom=DONUT_PANEL_MARGIN_BOTTOM_PX,
    )

    cns.barplot(data=df, y=category_column, x=percentage_column, ax=ax_bar)
    ax_bar.set_ylabel(ylabel)
    ax_bar.set_xlabel(xlabel)
    ax_bar.set_title(title)
    ax_bar.set_xlim(0, 100)
    for index, row in df.reset_index(drop=True).iterrows():
        ax_bar.text(
            row[percentage_column], index, f"{row[percentage_column]:.1f}%",
            ha="left", va="center",
        )

    # Every mp.panel() call is done, so the host's rectangle is final. Draw
    # first: before the first draw the panel reserves are unmeasured and
    # get_position() returns a provisional value.
    fig = multipanel.fig
    fig.canvas.draw()
    host_box = host.get_position()
    host.set_axis_off()
    host.set_gid(DONUT_GRID_HOST_GID)

    gridspec = fig.add_gridspec(
        n_donut_rows,
        DONUTS_PER_ROW,
        left=host_box.x0,
        right=host_box.x1,
        bottom=host_box.y0,
        top=host_box.y1,
        wspace=DONUT_GRID_WSPACE,
        hspace=DONUT_GRID_HSPACE,
    )

    # Donut cells: part vs whole per category, in the frame's row order.
    # All cells get their own letter from labels[1] onwards.
    for cell_index, (label, (_, row)) in enumerate(zip(labels[1:], df.iterrows(), strict=True)):
        ax_donut = fig.add_subplot(
            gridspec[cell_index // DONUTS_PER_ROW, cell_index % DONUTS_PER_ROW]
        )
        cns.setup_ax(ax_donut)
        plt.sca(ax_donut)
        cns.add_panel_label(label)

        category = str(row[category_column])
        part, whole = int(row[part_column]), int(row[whole_column])
        if part + whole == 0:
            ax_donut.set_axis_off()
            ax_donut.text(0.5, 0.5, f"{category}\n(no data)", ha="center", va="center",
                          transform=ax_donut.transAxes)
            continue

        draw_donut_panel(
            ax_donut, part, whole,
            part_label=part_label, whole_label=whole_label,
            center_text=f"{_wrap_center_text(category)}\n(n = {part + whole:,})",
        )

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
