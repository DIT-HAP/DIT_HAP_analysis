#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gene-Group Coherence — Visualization
====================================

Renders the overview figure for compute_coherence.py: the group-size and
z-score distributions, the centroid map in normalized fitness space, and — when
the metrics table carries them — the coherence-versus-biology panels.

This is the plotting companion, split out per ADR-0001. It computes nothing:
the shared-subunit fraction and the abundance/conservation uniformity terms are
columns of the metrics Parquet, written by compute_coherence.py, so the figure is
a pure renderer and every number it draws can be read back from the table.

Input
-----
- coherence_metrics.parquet: per-group metrics. Only these columns are read:
  n_scored_members, median_pairwise_distance_z, geom_median_DR, geom_median_DL,
  group_name, q_value, and whichever of frac_shared_members / abundance_cv /
  conservation_cv are present. Those last three gate the biology panels: an
  absent (or all-NaN) column drops its panel rather than drawing an empty one.

Output
------
- coherence.pdf (+ a .review.png sibling via save_dual)

Usage
-----
    python plot_coherence.py \\
        --input results/3a_coherence/{dataset}/{source}/coherence_metrics.parquet \\
        --output results/3a_coherence/{dataset}/{source}/coherence.pdf

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-09-03
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any

# 2. Data Processing Imports
import cnsplots as cns
import numpy as np
import pandas as pd
from matplotlib.axes import Axes

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from figure_render.histogram import draw_histogram_panel  # noqa: E402
from figures import apply_house_style, house_colors, save_dual  # noqa: E402
from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
_SQUARE_WIDTH, _SQUARE_HEIGHT = (100, 100)  # PanelShape.SQUARE, in layout pixels

# The centroid panel keeps the SQUARE footprint like its neighbours; its size
# legend and colourbar live OUTSIDE in the reserved right margin instead of
# widening the panel, so every panel in the row is the same size.
_CENTROID_MARGIN_RIGHT = 74

# Panel footprint = axes box + the panel label's pad + the inter-panel margin.
# multipanel wraps a row the moment its running width would exceed max_width, and
# does so silently, so the page is derived from the widest row instead of guessed.
_PANEL_MARGIN = 10

# Vertical space reserved under a panel whose row has another row beneath it.
# multipanel refines only the left and top decorations after a draw, so a bottom
# row's xlabel and tick labels are not measured and will overlap the next row's
# title unless reserved here.
_ROW_GAP = 46

# Inset colourbar geometry as an axes-fraction (x0, y0, width, height). x0 > 1
# puts it in the panel's reserved right margin, i.e. outside the plot area, below
# the size legend that take_legend_out anchors at the margin's top.
#
# x0 is what centres the colourbar under that legend. Matplotlib anchors the legend
# by a corner and this inset by the bar's left edge, so equal left offsets make the
# two blocks read as left-aligned even though they differ in width (the legend's
# title is wider than its entries; the bar carries tick labels and a rotated axis
# label). 1.147 shifts the bar so the two block centres coincide — measured at
# 17.3 px on a 200 px-wide axes, and identical across go_macrocomplex / go_cc /
# go_bp because every tick label is a single signed digit.
_CBAR_BOUNDS = (1.147, 0.02, 0.045, 0.42)

# Biology panels, in draw order. A 100 px panel fits about 25 characters at the
# house title size, so titles stay short — the y-axis already says "z-score", so
# the "coherence vs" half of the old titles was restating it and overflowing the
# panel into its neighbour.
# (column, x label, title)
_BIOLOGY_PANELS = [
    ("frac_shared_members", "Shared-subunit fraction", "Shared subunits"),
    ("abundance_cv", "Abundance CV", "Abundance uniformity"),
    ("conservation_cv", "Conservation CV", "Conservation uniformity"),
]

# Representative group sizes for the centroid map's size legend.
_SIZE_LEGEND_VALUES = [3, 10, 30, 100]
# Point area for a group of size s: log-compressed so a 300-member term does not
# swamp a 3-member one, scaled to be legible in a 150 px panel.
_SIZE_SCALE = 20.0

_REQUIRED_COLUMNS = ["n_scored_members", "median_pairwise_distance_z",
                     "geom_median_DR", "geom_median_DL", "group_name", "q_value"]

# Panels are lettered explicitly so the FDR panels' letters are predictable whatever
# number of biology columns the table carries: row 1 is A/B/C, the biology panels sit
# under them as D/E/F, and the FDR panels take the letters after those.
_BIOLOGY_LETTERS = ("D", "E", "F")

# cnsplots' own diverging scale, for the signed z-score. Resolved through
# cns.palettes() because it is not registered with matplotlib's cmap registry.
_DIVERGING_CMAP = "BuRd_custom"

# --- FDR-vs-coherence panel -------------------------------------------------
# x is -log10(q), not q: q spans 0.02..1 with a median near 0.1, so a linear axis
# crushes every interesting point against the left edge, and the significance
# boundary becomes a vertical line at -log10(q_max) instead of a judgement call.
# --- FDR-vs-coherence panels ------------------------------------------------
# Two panels over the SAME data with different x encodings, side by side, so the two
# axis choices can be compared before one is dropped.
#
# -log10(q) spreads a q that spans 0.02..1 with a median near 0.1; on a linear axis
# every interesting point is crushed against the left edge. Plain q is what a reader
# takes at face value. Whichever survives the comparison keeps its entry here and the
# other goes.
_FDR_PANELS = (
    ("_log10_q", "-log10(q)"),
    ("_q", "q (BH-adjusted)"),
)

# Floor for -log10(q). q is strictly positive by construction (the add-one
# estimator), but a pathologically tiny value would stretch the axis to nothing.
_FDR_Q_FLOOR = 1e-12

# The two FDR panels sit side by side and are sized so row 3 lines up with row 2:
# same left edge and same total width as the three SQUARE panels above. Those
# 3-column rows are what set the page's visual frame, and a row 3 that starts
# further right or runs wider reads as belonging to a different figure. Measured
# (go_macrocomplex, and identical for go_cc / go_bp): at width 159 with no left
# margin, G's axes left edge is 25.0 px — exactly D's — and the row spans 734.6 px
# against the SQUARE row's 735.6 px. The starting point was 200 px wide with a 14 px
# left margin, which put G 28 px right of D and made the row 191 px too wide.
#
# Width is a measured constant rather than a derived one because a panel's footprint
# includes the label and axis reserves multipanel only knows after it has drawn.
#
# The left margin is 0 on purpose: multipanel starts every row at the same x and each
# panel then adds its own margin_left, so a nonzero one is exactly what pushed row 3
# to the right. The right margin only sets the G-to-H gap.
#
# Height equals width: a 159 px-tall panel holds the name-label stacks without
# crowding (tightest inter-label gap 3.8 px, measured), and square matches the SQUARE
# panels' own aspect instead of making row 3 the one row of letterbox panels.
_FDR_PANEL_WIDTH, _FDR_PANEL_HEIGHT = 159, 159
_FDR_MARGIN_RIGHT = 14

# Labels sit INSIDE the axes, in the two empty quadrants of the S curve, with a leader
# line back to the point. Both columns meet at the same mid-panel x: the incoherent
# points (upper left) take their text rightwards into the upper-right gap, the coherent
# points (lower right) take theirs leftwards into the lower-left one. Placing them
# outside the frame instead would need a wide margin and drag every leader across the
# data.
#
# Each column is confined to its own quadrant so the two never share a y band, and the
# stack gap is computed from the wrapped line count — a fixed gap overlaps as soon as a
# name needs three lines.
#
# The two vertical metrics are DERIVED from the panel height rather than fixed axes
# fractions, because text is sized in points while these positions are axes fractions:
# the same 5 pt label eats a bigger slice of a shorter panel, so a constant tuned at one
# height silently crowds at another. 1 layout px is 1 pt (multipanel sizes the figure as
# px / 72 inches), so a line of text is _LABEL_FONT_SIZE / height in axes fraction and
# the multipliers below are just line spacing and clearance in units of the font size
# (1.6 and 1.2 reproduce the values that were hand-tuned at height 240).
_LABEL_COLUMN = 0.47
_LABEL_FONT_SIZE = 5
_LABEL_WRAP_WIDTH = 26        # characters per line before wrapping
_LABEL_LINE_HEIGHT = 1.6 * _LABEL_FONT_SIZE / _FDR_PANEL_HEIGHT   # one rendered line
_LABEL_BLOCK_GAP = 1.2 * _LABEL_FONT_SIZE / _FDR_PANEL_HEIGHT     # between two blocks
_LABEL_PADDING = 0.05         # how far past the extreme point the stack may reach
_LABEL_QUADRANTS = {"right": (0.52, 0.98), "left": (0.02, 0.48)}

# LabelSettings defaults, shared by the dataclass fields and argparse. They cannot
# be read off the dataclass: with slots=True, `LabelSettings.q_max` is a
# member_descriptor rather than the default value.
DEFAULT_LABEL_Q_MAX = 0.05
DEFAULT_LABEL_QUANTILE = 0.05
DEFAULT_LABEL_MAX = 5

# Page width, in layout pixels. 500 is measured, not derived: the widest row (the two
# FDR panels plus their margins) needs 481, and a fourth SQUARE panel on row 1 would
# need 568. Anything in between keeps the 3 x 2 layout intact; anything larger and
# multipanel silently packs that fourth panel onto row 1.
_MAX_WIDTH = 500


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class LabelSettings:
    """Which groups the FDR panel names: a significance cutoff and a per-side cap."""
    q_max: float = DEFAULT_LABEL_Q_MAX     # q at or below this counts as significant
    quantile: float = DEFAULT_LABEL_QUANTILE   # per-side share of the significant groups to name
    max_labels: int = DEFAULT_LABEL_MAX    # hard cap per side, so a 1,400-group source cannot flood the margin


@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Inputs and outputs for coherence visualization."""
    input_metrics: Path
    output: Path
    labels: LabelSettings = LabelSettings()

    def validate(self) -> None:
        """Raise ValueError if inputs are missing, then create output dirs."""
        if not self.input_metrics.exists():
            raise ValueError(f"Required input not found: {self.input_metrics}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
def biology_panels(table: pd.DataFrame) -> list[tuple[str, str, str]]:
    """The biology panels this table can support, in draw order."""
    # A column qualifies only if it is present and carries at least one value:
    # compute_coherence.py omits the feature CVs entirely when no features table
    # was passed, and leaves NaN for groups whose members lacked the feature, so
    # "absent" and "all-NaN" both mean there is nothing to draw.
    return [
        (column, xlabel, title)
        for column, xlabel, title in _BIOLOGY_PANELS
        if column in table.columns and table[column].notna().any()
    ]


def point_sizes(n_scored_members: pd.Series) -> np.ndarray:
    """Marker area per group from its scored-member count (log-compressed)."""
    return np.log1p(n_scored_members.to_numpy(dtype=float)) * _SIZE_SCALE


# =============================================================================
# CORE LOGIC — FDR panel labelling
# =============================================================================
def spread_positions(values: list[float], low: float, high: float, gap: float) -> list[float]:
    """Evenly de-overlap ascending positions while keeping them inside [low, high]."""
    # The classic label-spreading pass: a forward sweep pushes each label clear of
    # the one above it, then a backward sweep pulls the stack back inside the axis.
    # When the labels cannot all fit, the gap is relaxed rather than letting them
    # march off the top.
    count = len(values)
    if count == 0:
        return []
    span = high - low
    gap = min(gap, span / (count - 1)) if count > 1 else 0.0

    spread = [min(max(value, low), high) for value in values]
    for index in range(1, count):
        spread[index] = max(spread[index], spread[index - 1] + gap)
    spread[-1] = min(spread[-1], high)
    for index in range(count - 2, -1, -1):
        spread[index] = min(spread[index], spread[index + 1] - gap)
    # The sweeps accumulate float error, so a clamped end can land an ulp outside.
    return [min(max(value, low), high) for value in spread]


def label_count(available: int, quantile: float, max_labels: int) -> int:
    """How many groups to name on one side: the quantile share, capped at max_labels."""
    return min(max_labels, int(np.ceil(quantile * available)))


def labelled_extremes(
    table: pd.DataFrame, side: str, quantile: float, q_max: float, max_labels: int
) -> pd.DataFrame:
    """The groups to name at one end of the z-score axis, most extreme first."""
    # `coherent` takes the most negative z among the FDR-significant groups.
    # `incoherent` CANNOT use FDR: the coherence p-value is one-sided for tightness,
    # so a group more dispersed than random has p near 1 by construction (measured:
    # every q<=0.05 group in all three sources has z < 0). Selecting that end by
    # significance would return nothing, so it is selected by z alone.
    ordered = table.sort_values("median_pairwise_distance_z")
    if side == "coherent":
        eligible = ordered[ordered["q_value"] <= q_max]
        return eligible.head(label_count(len(eligible), quantile, max_labels))
    incoherent = ordered[ordered["median_pairwise_distance_z"] > 0]
    return incoherent.tail(label_count(len(incoherent), quantile, max_labels)).iloc[::-1]


def label_block_gap(rows: pd.DataFrame) -> float:
    """Vertical gap one row of labels needs, from how many lines its names wrap to."""
    lines = max(
        len(textwrap.wrap(str(name), width=_LABEL_WRAP_WIDTH)) for name in rows["group_name"]
    )
    return lines * _LABEL_LINE_HEIGHT + _LABEL_BLOCK_GAP


def annotate_extremes(ax: Axes, rows: pd.DataFrame, direction: str, x_column: str) -> None:
    """Name each row with a leader line into the panel's empty quadrant."""
    # Labels live INSIDE the axes: the S curve leaves the upper-right and lower-left
    # quadrants empty, so each column is confined to its own quadrant and the two never
    # share a y band. Rows and stack positions are both ordered by z ascending, so the
    # leaders fan out without crossing.
    #
    # The gap comes from the wrapped line count rather than a constant: a fixed gap
    # overlaps the moment a name needs a third line.
    if rows.empty:
        return
    rows = rows.sort_values("median_pairwise_distance_z")
    y_limits = ax.get_ylim()
    fractions = [(value - y_limits[0]) / (y_limits[1] - y_limits[0])
                 for value in rows["median_pairwise_distance_z"]]

    quadrant_low, quadrant_high = _LABEL_QUADRANTS[direction]
    gap = label_block_gap(rows)
    low = max(quadrant_low, fractions[0] - _LABEL_PADDING)
    high = min(quadrant_high, fractions[-1] + _LABEL_PADDING)
    if high - low < gap * (len(fractions) - 1):
        # Too many labels for the band their points occupy; use the whole quadrant and
        # let spread_positions relax the gap further if even that is not enough.
        low, high = quadrant_low, quadrant_high

    positions = spread_positions(fractions, low, high, gap)
    for (_, row), y_fraction in zip(rows.iterrows(), positions):
        ax.annotate(
            # fill() wraps rather than shortens: a reader needs the whole term name,
            # and the quadrant has room for two or three lines.
            textwrap.fill(str(row["group_name"]), width=_LABEL_WRAP_WIDTH),
            xy=(row[x_column], row["median_pairwise_distance_z"]),
            xytext=(_LABEL_COLUMN, y_fraction), textcoords="axes fraction",
            ha="left" if direction == "right" else "right", va="center",
            fontsize=_LABEL_FONT_SIZE,
            arrowprops={"arrowstyle": "-", "color": cns.GRAY, "linewidth": 0.6,
                        "shrinkA": 2, "shrinkB": 3},
        )


def plot_coherence(table: pd.DataFrame, settings: LabelSettings | None = None) -> None:
    """Draw the coherence overview onto a fresh house-styled multipanel figure."""
    # Panels, in order: group-size distribution, z-score distribution, centroid
    # map, then one per available biology column. multipanel labels every panel
    # automatically, A onwards — which is what fixes the old figure's A, B, D
    # (a panel C was dropped when the STRING channel was removed, and the letters
    # were never renumbered).
    apply_house_style()
    settings = settings or LabelSettings()

    if table.empty:
        ax = cns.multipanel().panel(width=_SQUARE_WIDTH, height=_SQUARE_HEIGHT)
        ax.text(0.5, 0.5, "No groups passed the size filter", ha="center", va="center")
        ax.set_axis_off()
        return

    # multipanel sizes each panel from its own rendered decorations, so panels
    # whose labels are filled in are created with an explicit pad rather than
    # relying on a grid to align them.
    #
    # Rows wrap on width, so max_width is what holds the layout at 3 panels per row:
    # row 1 (A, B, C plus C's legend margin) is 413 layout px, and a fourth SQUARE
    # panel would need ~122 more, so 500 wraps it. Raise max_width past ~535 and
    # multipanel silently packs four panels onto row 1 and the grid becomes 4 x 2.
    multipanel = cns.multipanel(max_width=_MAX_WIDTH)

    ax_size = multipanel.panel("A", width=_SQUARE_WIDTH, height=_SQUARE_HEIGHT, margin_bottom=_ROW_GAP)
    draw_histogram_panel(
        ax_size, table["n_scored_members"], bins=21, log_scale=True,
        xlabel="Group size\n(DR<threshold members)", ylabel="Number of groups",
        title="Group size distribution",
    )

    ax_z = multipanel.panel("B", width=_SQUARE_WIDTH, height=_SQUARE_HEIGHT, margin_bottom=_ROW_GAP)
    draw_histogram_panel(
        ax_z, table["median_pairwise_distance_z"], bins=20,
        xlabel="MPD z-score\n(negative = coherent)", ylabel="Number of groups",
        title="Coherence z-scores",
    )
    # z = 0 is the null: at or above it the group is no tighter than a random
    # draw of the same size. Drawn after the histogram so it reads on top.
    ax_z.axvline(0.0, color=cns.GRAY, linestyle="--", linewidth=1.0)

    ax_centroid = multipanel.panel(
        "C", width=_SQUARE_WIDTH, height=_SQUARE_HEIGHT,
        margin_right=_CENTROID_MARGIN_RIGHT, margin_bottom=_ROW_GAP,
    )
    sizes = point_sizes(table["n_scored_members"])
    # z-score is signed and 0 is its meaningful midpoint, so a diverging map is the
    # right encoding, with vmin/vmax set symmetrically or the midpoint colour lands
    # at an arbitrary value. cnsplots ships this one (BuRd_custom) as its diverging
    # map; it runs blue at the low end, which puts coherent (negative z) on the cool
    # end. A plain matplotlib name like coolwarm would work but is not the
    # package's own scale.
    #
    # The edge is not decoration: a diverging map spends its midpoint on white, and
    # most groups sit at z ~ 0, so without a hairline outline they vanish into the
    # page.
    z_limits = float(np.nanmax(np.abs(table["median_pairwise_distance_z"].to_numpy(dtype=float)))) or 1.0
    scatter = ax_centroid.scatter(
        table["geom_median_DR"], table["geom_median_DL"],
        c=table["median_pairwise_distance_z"], s=sizes,
        cmap=cns.palettes(_DIVERGING_CMAP), vmin=-z_limits, vmax=z_limits,
        alpha=0.85, edgecolors=cns.GRAY, linewidths=0.3,
    )
    ax_centroid.set(
        xlabel="typical DR", ylabel="typical DL/10", title="Group centroid positions"
    )

    colorbar_ax = ax_centroid.inset_axes(_CBAR_BOUNDS)
    colorbar = ax_centroid.figure.colorbar(scatter, cax=colorbar_ax)
    colorbar.ax.tick_params(length=0, pad=1)
    colorbar.set_label("z-score", labelpad=1)

    legend_handles = [
        ax_centroid.scatter([], [], s=float(np.log1p(value) * _SIZE_SCALE),
                            color=cns.GRAY, alpha=0.7, edgecolors="none")
        for value in _SIZE_LEGEND_VALUES
    ]
    ax_centroid.legend(
        legend_handles, [str(value) for value in _SIZE_LEGEND_VALUES], frameon=False,
    )
    # Keyword-only: written positionally this silently sets the legend TITLE.
    cns.take_legend_out("Group size", ax=ax_centroid)

    # Every biology panel plots z-score on y, so they share one y range: with
    # autoscaled panels a group at z=-6 renders exactly like one at z=-1, and the
    # panels only mean anything compared against each other. Their x variables are
    # different quantities on different scales, so x stays per-panel.
    biology = biology_panels(table)
    z_values = table["median_pairwise_distance_z"].to_numpy(dtype=float)
    shared_z_limits = (float(np.nanmin(z_values)), float(np.nanmax(z_values)))

    for letter, (column, xlabel, title) in zip(_BIOLOGY_LETTERS, biology):
        ax = multipanel.panel(
            letter, width=_SQUARE_WIDTH, height=_SQUARE_HEIGHT, margin_bottom=_ROW_GAP
        )
        cns.scatterplot(table, column, "median_pairwise_distance_z", ax=ax,
                        color=house_colors((3,))[0])
        ax.set(xlabel=xlabel, ylabel="z-score", title=title)
        ax.set_ylim(*shared_z_limits)
        ax.axhline(0.0, color=cns.GRAY, linestyle="--", linewidth=1.0)

    # The FDR row follows whichever biology panels were drawn, so its letters depend
    # on how many there are.
    draw_fdr_panels(multipanel, table, settings, first_letter=chr(ord(_BIOLOGY_LETTERS[0]) + len(biology)))


def fdr_axis(q_values: pd.Series, encoding: str, q_max: float) -> tuple[pd.Series, float]:
    """X values for an FDR panel and the significance cutoff on the same scale."""
    if encoding == "_log10_q":
        return -np.log10(q_values.clip(lower=_FDR_Q_FLOOR)), -np.log10(q_max)
    return q_values, q_max


def draw_fdr_panels(
    multipanel: Any, table: pd.DataFrame, settings: LabelSettings, first_letter: str
) -> None:
    """Draw the FDR-versus-coherence scatter once per x encoding, side by side."""
    for offset, (encoding, xlabel) in enumerate(_FDR_PANELS):
        ax = multipanel.panel(
            chr(ord(first_letter) + offset),
            width=_FDR_PANEL_WIDTH, height=_FDR_PANEL_HEIGHT,
            margin_right=_FDR_MARGIN_RIGHT,
        )
        x_values, cutoff = fdr_axis(table["q_value"], encoding, settings.q_max)
        plotted = table.assign(**{encoding: x_values})
        cns.scatterplot(plotted, encoding, "median_pairwise_distance_z", ax=ax,
                        color=house_colors((3,))[0])
        ax.set(xlabel=xlabel, ylabel="z-score", title="Coherence vs significance")
        ax.axhline(0.0, color=cns.GRAY, linestyle="--", linewidth=1.0)
        ax.axvline(cutoff, color=cns.GRAY, linestyle=":", linewidth=1.0)

        # Each end is named into the quadrant its points leave empty. The incoherent
        # ones sit upper-left and the coherent ones lower-right, because the coherence
        # p-value is one-sided for tightness: no group more dispersed than random can
        # be FDR-significant.
        for end, direction in (("incoherent", "right"), ("coherent", "left")):
            rows = labelled_extremes(plotted, end, settings.quantile, settings.q_max, settings.max_labels)
            annotate_extremes(ax, rows, direction, encoding)


# =============================================================================
# MAIN EXECUTION
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotConfig) -> None:
    """Generate the coherence overview figure from the computed metrics."""
    config.validate()
    table = read_parquet(config.input_metrics)
    missing = [column for column in _REQUIRED_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(
            f"metrics table missing required column(s) {missing} (have: {list(table.columns)})"
        )

    logger.info(
        f"Loaded {len(table):,} groups; biology panels: "
        f"{[column for column, _, _ in biology_panels(table)] or 'none'}"
    )
    plot_coherence(table, config.labels)
    save_dual(config.output.with_suffix(""))
    logger.success(f"Wrote {config.output}")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Generate coherence visualization")
    parser.add_argument("--input", dest="input_metrics", type=Path, required=True,
                        help="Input coherence metrics Parquet from compute_coherence.py")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output coherence figure PDF")
    parser.add_argument("--label-q-max", type=float, default=DEFAULT_LABEL_Q_MAX,
                        help="FDR panel: q at or below this counts as significant")
    parser.add_argument("--label-quantile", type=float, default=DEFAULT_LABEL_QUANTILE,
                        help="FDR panel: per-side share of significant groups to name")
    parser.add_argument("--label-max", type=int, default=DEFAULT_LABEL_MAX,
                        help="FDR panel: hard cap on names per side")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run plotting, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotConfig(
            input_metrics=args.input_metrics,
            output=args.output,
            labels=LabelSettings(
                q_max=args.label_q_max,
                quantile=args.label_quantile,
                max_labels=args.label_max,
            ),
        )
        run(config)
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
