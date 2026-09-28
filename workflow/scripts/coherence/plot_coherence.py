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
- coherence_metrics.parquet (or any table with the same columns): per-group
  metrics. Only these columns are read: n_scored_members,
  median_pairwise_distance_z, geom_median_DR, geom_median_DL, group_name,
  q_value, and whichever of frac_shared_members / abundance_cv /
  conservation_cv are present. Those last three gate the biology panels: an
  absent (or all-NaN) column drops its panel rather than drawing an empty one.

Two optional modes change only the colour encoding, never the panel layout, so a
comparison figure can be read panel-for-panel against a per-source one:

  --color-by source   Draw every panel once per `source` and colour by it. Panel C
                      gives up its z-score map for this (z is already the y axis of
                      every other panel); panels A/B become overlaid step outlines
                      on shared bin edges; one legend on panel A keys the figure.
  --dedup-series PATH Append a second table as an extra `dedup` series — the
                      de-duplicated representative subset, drawn alongside the full
                      sets rather than instead of them.

Output
------
- coherence.pdf (+ a .review.png sibling via save_dual)

Usage
-----
    python plot_coherence.py \\
        --input results/3a_coherence/{dataset}/{source}/coherence_metrics.parquet \\
        --output results/3a_coherence/{dataset}/{source}/coherence.pdf

    # The de-duplicated representative set on its own.
    python plot_coherence.py \\
        --input results/3a_coherence/{dataset}/coherence_terms_representatives.tsv \\
        --output results/3a_coherence/{dataset}/coherence_dedup.pdf

    # Every source + the representatives, coloured by source.
    python plot_coherence.py \\
        --input results/3a_coherence/{dataset}/coherence_metrics_combined.parquet \\
        --color-by source \\
        --dedup-series results/3a_coherence/{dataset}/coherence_terms_representatives.tsv \\
        --output results/3a_coherence/{dataset}/coherence_by_source.pdf

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-09-03
Version:  2.1.0
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
from figures import apply_house_style, apply_log_scale, house_colors, save_dual  # noqa: E402
from coherence.palette import SOURCE_ORDER, source_colors  # noqa: E402
from io_table import read_file  # noqa: E402
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

# --- Cross-source comparison mode (color_by="source") -----------------------
# Bin counts for panels A and B. Same as the single-series figure's, so the two can
# be read against each other; the comparison lays the per-source histograms as step
# OUTLINES on SHARED edges, because binning each source on its own range would make
# the panels compare bin widths as much as counts, and six overlapping filled bars
# would hide one another.
_COMPARISON_SIZE_BINS = 21
_COMPARISON_Z_BINS = 20

_REQUIRED_COLUMNS = ["n_scored_members", "median_pairwise_distance_z",
                     "geom_median_DR", "geom_median_DL", "group_name", "q_value"]

# Panels are lettered explicitly so the FDR panels' letters are predictable whatever
# number of biology columns the table carries: row 1 is A/B/C, the biology panels sit
# under them as D/E/F, and the FDR panels take the letters after those.
_BIOLOGY_LETTERS = ("D", "E", "F")

# cnsplots' own diverging scale, for the signed z-score. Resolved through
# cns.palettes() because it is not registered with matplotlib's cmap registry.
_DIVERGING_CMAP = "BuRd_custom"

# Half-range of the diverging colour scale, as a percentile of |z|.
#
# Scaling the map to max|z| is what made the panel look washed out: the z
# distribution is strongly asymmetric (median -1.5, |z| p98 ~3.2, but the tail
# runs to 6), so a symmetric +/-6 range spends most of its span on a handful of
# points and paints the bulk of the groups within half a step of the white
# midpoint. Clipping at a high percentile keeps the colour range where the data
# actually is; the colourbar's end triangles still say the tails were clipped.
_PANEL_C_PERCENTILE = 98.0

# Blended alpha for the centroid points. Below 1 the white midpoint of the
# diverging map lets the page through and washes out again, which is the second
# half of the same problem; the grey hairline outline still separates overlaps.
_PANEL_C_ALPHA = 1.0

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

# Labels sit INSIDE the axes, anchored beside their OWN point (a short leader),
# and are placed by scoring a grid of candidate anchors and keeping the cheapest.
#
# No single fixed rule places them: the same end's points sit at the left edge in
# one panel and the right in the other (q and -log10(q) order the x axis in
# opposite directions), so which side has room, and how far a label must travel to
# clear the curve, flips between the two. A shared mid-panel column — the previous
# approach — made every leader as long as its point sat far from the column, and in
# panel G those leaders ran the width of the panel and crossed the data.
#
# Scoring counts drawn points and already-placed labels inside a candidate's
# estimated footprint, so a label lands in white space with a leader no longer than
# it has to be. The y bands stay disjoint between the two ends, so they cannot
# collide with each other whatever the search picks.
#
# The two vertical metrics are DERIVED from the panel height rather than fixed axes
# fractions, because text is sized in points while these positions are axes fractions:
# the same 5 pt label eats a bigger slice of a shorter panel, so a constant tuned at one
# height silently crowds at another. 1 layout px is 1 pt (multipanel sizes the figure as
# px / 72 inches), so a line of text is _LABEL_FONT_SIZE / height in axes fraction and
# the multipliers below are just line spacing and clearance in units of the font size
# (1.6 and 1.2 reproduce the values that were hand-tuned at height 240).
_LABEL_FONT_SIZE = 5
_LABEL_WRAP_WIDTH = 35        # characters per line before wrapping
_LABEL_LINE_HEIGHT = 1.6 * _LABEL_FONT_SIZE / _FDR_PANEL_HEIGHT   # one rendered line
_LABEL_BLOCK_GAP = 1.2 * _LABEL_FONT_SIZE / _FDR_PANEL_HEIGHT     # between two blocks
_LABEL_ANCHOR_PAD = 0.035     # smallest gap between a point and its own label

# The y band each end's labels occupy. The incoherent end (z > 0) sits in the
# upper band and the coherent end (z < 0) in the lower one, because the coherence
# p-value is one-sided for tightness: no group more dispersed than random can be
# FDR-significant, so the S curve leaves the upper-right and lower-left quadrants
# empty. Disjoint bands are what keep the two stacks apart.
_LABEL_BANDS = {"incoherent": (0.52, 0.98), "coherent": (0.02, 0.48)}

# Placement search. Character width at _LABEL_FONT_SIZE, in axes fraction of a
# _FDR_PANEL_WIDTH-wide panel (~0.5 em per glyph for the house font) — used only to
# estimate a candidate's footprint for the overlap test; matplotlib lays out the
# text itself.
_LABEL_CHAR_WIDTH = 0.017
_POINT_RADIUS = 0.012         # drawn marker radius, same units
_ANCHOR_PUSHES = (0.0, 0.08, 0.16, 0.24, 0.34)   # how far out from the point to try
_ANCHOR_Y_STEPS = 24          # candidate rows within the end's band

# Scoring, in strict priority order: a candidate must never overlap another label's
# text, may cover a few markers rather than none, and only then is a leader that cuts
# across a label penalised. Each later term's worst case is kept below the one before
# it, so a candidate can never buy its way out of the harder defect with a cheaper
# one — and the point term saturates, so a candidate in a dense region is not scored
# as impossible and the search does not abandon a clear row to escape it.
_POINT_HIT_PENALTY = 20.0
_POINT_HIT_CAP = 4            # x the penalty above = 80, below the label-overlap cost
_LEADER_CROSS_PENALTY = 20.0  # x four other labels = 80, likewise
_LABEL_HIT_PENALTY = 200.0    # text on text: the one outcome to avoid outright
_OVERLAP_AREA_WEIGHT = 5.0    # tie-break between candidates that all collide
_LEADER_PENALTY = 1.0         # per unit of leader length

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
    # "none" (default) = today's single-series figure with z encoded as colour.
    # "source" = the cross-source comparison: every panel is drawn once per `source`
    # and coloured by it, so the per-source figures can be read side by side.
    color_by: str = "none"
    # Appended as an extra `source` ("dedup") when set — the de-duplicated
    # representative subset, which is a subset of the rows already in
    # `input_metrics` and is drawn alongside them rather than instead of them.
    dedup_series: Path | None = None

    def validate(self) -> None:
        """Raise ValueError if inputs are missing, then create output dirs."""
        if not self.input_metrics.exists():
            raise ValueError(f"Required input not found: {self.input_metrics}")
        if self.dedup_series is not None and not self.dedup_series.exists():
            raise ValueError(f"Required input not found: {self.dedup_series}")
        if self.color_by not in ("none", "source"):
            raise ValueError(f"color_by must be 'none' or 'source': {self.color_by!r}")
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


def visible_fractions(
    ax: Axes, table: pd.DataFrame, x_column: str
) -> list[tuple[float, float]]:
    """Every drawn point as an axes-fraction (x, y) pair, for the placement search."""
    x_limits, y_limits = ax.get_xlim(), ax.get_ylim()
    x_span = x_limits[1] - x_limits[0]
    y_span = y_limits[1] - y_limits[0]
    return [
        ((x - x_limits[0]) / x_span, (z - y_limits[0]) / y_span)
        for x, z in zip(table[x_column], table["median_pairwise_distance_z"])
    ]


def boxes_hit(
    box: tuple[float, float, float, float], points: list[tuple[float, float]],
) -> int:
    """How many drawn points fall inside (or under) an axes-fraction box."""
    x0, x1, y0, y1 = box
    return sum(
        x0 - _POINT_RADIUS <= px <= x1 + _POINT_RADIUS
        and y0 - _POINT_RADIUS <= py <= y1 + _POINT_RADIUS
        for px, py in points
    )


def place_label(
    point: tuple[float, float],
    width: float,
    height: float,
    band: tuple[float, float],
    occupied: list[tuple[float, float, float, float]],
    points: list[tuple[float, float]],
) -> tuple[float, float, bool]:
    """Anchor and side for one label: the cheapest of a small candidate grid."""
    # Candidates differ in how far the label is pushed from its point (which picks
    # the side, by room) and in where it sits within the end's band. Each is scored
    # by the drawn points and placed labels it covers, the leaders it cuts, and the
    # leader length it costs — in that order of priority.
    #
    # Returns (anchor_x, anchor_y, to_the_left) so the caller can place the text and
    # record the box without re-deriving the side, which a pushed anchor would get
    # wrong near the middle of the panel.
    point_x, point_y = point
    to_the_left = point_x > 0.5
    low = band[0] + height / 2
    high = band[1] - height / 2
    if high < low:
        low = high = (band[0] + band[1]) / 2

    # The anchor is clamped rather than rejected, so a label wider than the room on
    # its side still gets a box inside the axes instead of falling back to the bare
    # point and sitting on top of the data.
    def anchor_at(push: float) -> float:
        offset = (push + _LABEL_ANCHOR_PAD) * (-1 if to_the_left else 1)
        raw = point_x + offset
        return min(max(raw, width), 1.0) if to_the_left else max(min(raw, 1.0 - width), 0.0)

    best, best_score = (point_x, point_y, to_the_left), None
    for push in _ANCHOR_PUSHES:
        anchor_x = anchor_at(push)
        for step in range(_ANCHOR_Y_STEPS + 1):
            anchor_y = low + (high - low) * step / _ANCHOR_Y_STEPS
            box = box_for(anchor_x, anchor_y, width, height, to_the_left)
            hits = boxes_hit(box, points)
            score = _LEADER_PENALTY * (abs(anchor_x - point_x) + abs(anchor_y - point_y))
            score += _POINT_HIT_PENALTY * min(hits, _POINT_HIT_CAP)
            score += _LABEL_HIT_PENALTY * sum(overlap_area(box, other) > 0.0 for other in occupied)
            score += _OVERLAP_AREA_WEIGHT * sum(overlap_area(box, other) for other in occupied)
            # A leader has to reach its label without cutting through another one, or
            # the panel reads as a tangle even when every box is clear of the rest.
            score += _LEADER_CROSS_PENALTY * leader_crossings((point_x, point_y), (anchor_x, anchor_y), occupied)
            if best_score is None or score < best_score:
                best, best_score = (anchor_x, anchor_y, to_the_left), score

    return best


def box_for(
    anchor_x: float, anchor_y: float, width: float, height: float, to_the_left: bool,
) -> tuple[float, float, float, float]:
    """The (x0, x1, y0, y1) axes-fraction footprint of a label anchored at (x, y)."""
    left, right = (anchor_x - width, anchor_x) if to_the_left else (anchor_x, anchor_x + width)
    return (left, right, anchor_y - height / 2, anchor_y + height / 2)


def overlap_area(
    box: tuple[float, float, float, float], other: tuple[float, float, float, float],
) -> float:
    """Intersection area of two boxes, 0.0 when they are clear of each other."""
    dx = min(box[1], other[1]) - max(box[0], other[0])
    dy = min(box[3], other[3]) - max(box[2], other[2])
    return max(dx, 0.0) * max(dy, 0.0)


def leader_crossings(
    point: tuple[float, float],
    anchor: tuple[float, float],
    occupied: list[tuple[float, float, float, float]],
) -> int:
    """How many placed labels a leader from `point` to `anchor` cuts through."""
    # Sampled rather than solved: 25 points along a segment a tenth of the panel
    # long cannot miss a label box, and the test stays a few lines.
    steps = 24
    hits = 0
    for box in occupied:
        for step in range(steps + 1):
            fraction = step / steps
            x = point[0] + (anchor[0] - point[0]) * fraction
            y = point[1] + (anchor[1] - point[1]) * fraction
            if box[0] <= x <= box[1] and box[2] <= y <= box[3]:
                hits += 1
                break
    return hits


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


def label_block_lines(rows: pd.DataFrame) -> int:
    """Rendered line count of the tallest name in `rows`."""
    return max(
        len(textwrap.wrap(str(name), width=_LABEL_WRAP_WIDTH)) for name in rows["group_name"]
    )


def annotate_extremes(
    ax: Axes, rows: pd.DataFrame, end: str, x_column: str,
    points: list[tuple[float, float]],
) -> None:
    """Name each row with a leader line out of its own point."""
    if rows.empty:
        return
    rows = rows.sort_values("median_pairwise_distance_z")
    x_limits, y_limits = ax.get_xlim(), ax.get_ylim()
    x_span, y_span = x_limits[1] - x_limits[0], y_limits[1] - y_limits[0]
    band = _LABEL_BANDS[end]
    lines = label_block_lines(rows)
    # The block gap is part of the footprint the search reserves, so two labels on
    # adjacent candidates come out spaced rather than flush.
    height = lines * _LABEL_LINE_HEIGHT + _LABEL_BLOCK_GAP

    entries = []
    for _, row in rows.iterrows():
        text = textwrap.fill(str(row["group_name"]), width=_LABEL_WRAP_WIDTH)
        entries.append((
            text,
            max(len(line) for line in text.splitlines()) * _LABEL_CHAR_WIDTH,
            ((row[x_column] - x_limits[0]) / x_span,
             (row["median_pairwise_distance_z"] - y_limits[0]) / y_span),
            (row[x_column], row["median_pairwise_distance_z"]),
        ))

    placements = place_all(entries, height, [band] * len(entries), points)
    if any_overlap(placements, entries, height):
        # The bands are only about four text blocks tall, so the search's greedy
        # picks can crowd out the labels still to come and leave two on top of each
        # other. An even stack over the band cannot overlap at all, so it takes over
        # rather than emitting a collision.
        rows_y = spread_positions(
            [point[1] for _, _, point, _ in entries],
            band[0] + height / 2, band[1] - height / 2, height,
        )
        placements = place_all(
            entries, height, [(y - height / 2, y + height / 2) for y in rows_y], points
        )

    for (text, _, _, xy), (anchor_x, anchor_y, to_the_left) in zip(entries, placements):
        ax.annotate(
            text,
            xy=xy, xytext=(anchor_x, anchor_y), textcoords="axes fraction",
            ha="right" if to_the_left else "left", va="center",
            fontsize=_LABEL_FONT_SIZE,
            arrowprops={"arrowstyle": "-", "color": cns.GRAY, "linewidth": 0.6,
                        "shrinkA": 2, "shrinkB": 3},
        )


def place_all(
    entries: list[tuple[str, float, tuple[float, float], tuple[float, float]]],
    height: float,
    bands: list[tuple[float, float]],
    points: list[tuple[float, float]],
) -> list[tuple[float, float, bool]]:
    """Run the placement search over every label, later ones avoiding earlier boxes."""
    occupied: list[tuple[float, float, float, float]] = []
    placements = []
    for (_, width, point, _), band in zip(entries, bands):
        anchor = place_label(point, width, height, band, occupied, points)
        occupied.append(box_for(anchor[0], anchor[1], width, height, anchor[2]))
        placements.append(anchor)
    return placements


def any_overlap(
    placements: list[tuple[float, float, bool]],
    entries: list[tuple[str, float, tuple[float, float], tuple[float, float]]],
    height: float,
) -> bool:
    """Whether the placed label boxes collide with each other anywhere."""
    boxes = [
        box_for(anchor[0], anchor[1], width, height, anchor[2])
        for anchor, (_, width, _, _) in zip(placements, entries)
    ]
    return any(
        overlap_area(a, b) > 0.0
        for index, a in enumerate(boxes)
        for b in boxes[index + 1:]
    )


def _source_series(table: pd.DataFrame) -> list[str]:
    """The sources present in `table`, in SOURCE_ORDER (anything unlisted last)."""
    # Fixed order, not order of appearance: the legend, the draw order and the
    # colour are then the same however the table happens to be sorted, and a source
    # keeps its colour between the two figures that use this module.
    present = set(table["source"].dropna().unique())
    ordered = [source for source in SOURCE_ORDER if source in present]
    return ordered + sorted(present.difference(ordered))


def draw_source_histogram(
    ax: Axes, table: pd.DataFrame, column: str, *, bins: int, log_scale: bool,
    colors: dict[str, str], series: list[str], legend: bool = False,
) -> None:
    """One step-histogram per source on SHARED bin edges."""
    # Shared edges come from the whole table, not per source: binning each source on
    # its own range would make two panels compare bin widths rather than counts.
    # Outlines, not the house's filled bars, because six overlapping fills are
    # unreadable; the single-series figure keeps its filled bar in panel A.
    values = table[[column, "source"]].dropna()
    if log_scale:
        values = values[values[column] > 0]
    if values.empty:
        logger.warning(f"Panel {column!r} has no valid data")
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
        return

    if log_scale:
        edges: np.ndarray | int = np.logspace(
            np.log10(values[column].min()), np.log10(values[column].max()), bins + 1
        )
    else:
        edges = bins

    for source in series:
        source_values = values.loc[values["source"] == source, column]
        if source_values.empty:
            continue
        ax.hist(source_values, bins=edges, histtype="step", linewidth=1.0,
                color=colors[source], label=source)

    if log_scale:
        # The edges are already laid out in log10 space, so the axis is switched
        # afterwards rather than by a log flag that would re-bin them.
        apply_log_scale(ax, x=True, y=False)
    if legend:
        ax.legend(frameon=False, fontsize=5, loc="upper right")


def draw_source_scatter(
    ax: Axes, table: pd.DataFrame, x_column: str, y_column: str, *,
    by_source: bool, colors: dict[str, str], series: list[str],
) -> None:
    """Scatter `x_column`-vs-`y_column`, one series per source in comparison mode."""
    # Goes through cns.scatterplot once per source rather than ax.scatter, so the
    # marker size and style stay exactly what the single-series panels draw — only
    # the colour changes.
    if not by_source:
        cns.scatterplot(table, x_column, y_column, ax=ax, color=house_colors((3,))[0])
        return
    for source in series:
        rows = table[table["source"] == source]
        if rows.empty:
            continue
        cns.scatterplot(rows, x_column, y_column, ax=ax, color=colors[source])


def plot_coherence(
    table: pd.DataFrame, settings: LabelSettings | None = None, color_by: str = "none"
) -> None:
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

    # Comparison mode: one flat colour per source, replacing panel C's z-score map.
    # The z it gives up is not lost from the figure — every biology and FDR panel
    # below already plots z on y — while what the comparison is for (which source a
    # group came from) has nowhere else to live.
    by_source = color_by == "source"
    series = _source_series(table) if by_source else []
    colors = source_colors(series) if by_source else {}

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
    if by_source:
        # The one legend for the whole figure: the key is read once, at the top.
        draw_source_histogram(
            ax_size, table, "n_scored_members", bins=_COMPARISON_SIZE_BINS, log_scale=True,
            colors=colors, series=series, legend=True,
        )
        ax_size.set(
            xlabel="Group size\n(DR<threshold members)", ylabel="Number of groups",
            title="Group size distribution",
        )
    else:
        draw_histogram_panel(
            ax_size, table["n_scored_members"], bins=21, log_scale=True,
            xlabel="Group size\n(DR<threshold members)", ylabel="Number of groups",
            title="Group size distribution",
        )

    ax_z = multipanel.panel("B", width=_SQUARE_WIDTH, height=_SQUARE_HEIGHT, margin_bottom=_ROW_GAP)
    if by_source:
        draw_source_histogram(
            ax_z, table, "median_pairwise_distance_z", bins=_COMPARISON_Z_BINS,
            log_scale=False, colors=colors, series=series,
        )
        ax_z.set(
            xlabel="MPD z-score\n(negative = coherent)", ylabel="Number of groups",
            title="Coherence z-scores",
        )
    else:
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
    # The limits are a percentile, not max|z| — see _PANEL_C_PERCENTILE. The
    # colourbar gets `extend="both"` so the clipped tails are visible as such.
    #
    # The edge is not decoration: a diverging map spends its midpoint on white, and
    # most groups sit at z ~ 0, so without a hairline outline they vanish into the
    # page.
    if by_source:
        # One flat colour per source; the z-score map and its colourbar are dropped.
        # Point size still carries group size, so the two encodings stay the same
        # shape as the per-source figure's C panel — only colour changes meaning.
        for source in series:
            rows = table[table["source"] == source]
            ax_centroid.scatter(
                rows["geom_median_DR"], rows["geom_median_DL"],
                s=point_sizes(rows["n_scored_members"]), color=colors[source],
                alpha=_PANEL_C_ALPHA, edgecolors=cns.GRAY, linewidths=0.3,
            )
    else:
        z_values = table["median_pairwise_distance_z"].to_numpy(dtype=float)
        z_limits = float(np.nanpercentile(np.abs(z_values), _PANEL_C_PERCENTILE))
        if not z_limits:
            z_limits = 1.0
        scatter = ax_centroid.scatter(
            table["geom_median_DR"], table["geom_median_DL"],
            c=table["median_pairwise_distance_z"], s=sizes,
            cmap=cns.palettes(_DIVERGING_CMAP), vmin=-z_limits, vmax=z_limits,
            alpha=_PANEL_C_ALPHA, edgecolors=cns.GRAY, linewidths=0.3,
        )
    ax_centroid.set(
        xlabel="typical DR", ylabel="typical DL/10", title="Group centroid positions"
    )

    if not by_source:
        colorbar_ax = ax_centroid.inset_axes(_CBAR_BOUNDS)
        colorbar = ax_centroid.figure.colorbar(scatter, cax=colorbar_ax, extend="both")
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
        draw_source_scatter(
            ax, table, column, "median_pairwise_distance_z",
            by_source=by_source, colors=colors, series=series,
        )
        ax.set(xlabel=xlabel, ylabel="z-score", title=title)
        ax.set_ylim(*shared_z_limits)
        ax.axhline(0.0, color=cns.GRAY, linestyle="--", linewidth=1.0)

    # The FDR row follows whichever biology panels were drawn, so its letters depend
    # on how many there are.
    draw_fdr_panels(
        multipanel, table, settings,
        first_letter=chr(ord(_BIOLOGY_LETTERS[0]) + len(biology)),
        by_source=by_source, colors=colors, series=series,
    )


def fdr_axis(q_values: pd.Series, encoding: str, q_max: float) -> tuple[pd.Series, float]:
    """X values for an FDR panel and the significance cutoff on the same scale."""
    if encoding == "_log10_q":
        return -np.log10(q_values.clip(lower=_FDR_Q_FLOOR)), -np.log10(q_max)
    return q_values, q_max


def draw_fdr_panels(
    multipanel: Any, table: pd.DataFrame, settings: LabelSettings, first_letter: str,
    by_source: bool = False, colors: dict[str, str] | None = None,
    series: list[str] | None = None,
) -> None:
    """Draw the FDR-versus-coherence scatter once per x encoding, side by side."""
    colors, series = colors or {}, series or []
    for offset, (encoding, xlabel) in enumerate(_FDR_PANELS):
        ax = multipanel.panel(
            chr(ord(first_letter) + offset),
            width=_FDR_PANEL_WIDTH, height=_FDR_PANEL_HEIGHT,
            margin_right=_FDR_MARGIN_RIGHT,
        )
        x_values, cutoff = fdr_axis(table["q_value"], encoding, settings.q_max)
        plotted = table.assign(**{encoding: x_values})
        draw_source_scatter(
            ax, plotted, encoding, "median_pairwise_distance_z",
            by_source=by_source, colors=colors, series=series,
        )
        ax.set(xlabel=xlabel, ylabel="z-score", title="Coherence vs significance")
        ax.axhline(0.0, color=cns.GRAY, linestyle="--", linewidth=1.0)
        ax.axvline(cutoff, color=cns.GRAY, linestyle=":", linewidth=1.0)

        # Each end is named beside its own points: the incoherent ones sit in the
        # upper band and the coherent ones in the lower, because the coherence
        # p-value is one-sided for tightness — no group more dispersed than random
        # can be FDR-significant.
        points = visible_fractions(ax, plotted, encoding)
        for end in ("incoherent", "coherent"):
            rows = labelled_extremes(plotted, end, settings.quantile, settings.q_max, settings.max_labels)
            annotate_extremes(ax, rows, end, encoding, points)


# =============================================================================
# MAIN EXECUTION
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotConfig) -> None:
    """Generate the coherence overview figure from the computed metrics."""
    config.validate()
    # read_file dispatches on extension: the per-source tables are Parquet and the
    # de-duplicated representatives (a final human-facing artifact) are TSV.
    table = read_file(config.input_metrics)
    missing = [column for column in _REQUIRED_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(
            f"metrics table missing required column(s) {missing} (have: {list(table.columns)})"
        )
    if config.dedup_series is not None:
        # Appended, not swapped in: the comparison puts the representative subset
        # next to the full per-source sets so the de-duplication's effect is visible
        # as a difference between two series on the same axes.
        dedup = read_file(config.dedup_series)
        missing = [column for column in _REQUIRED_COLUMNS if column not in dedup.columns]
        if missing:
            raise ValueError(
                f"dedup table missing required column(s) {missing} (have: {list(dedup.columns)})"
            )
        table = pd.concat([table, dedup.assign(source="dedup")], ignore_index=True)
    if config.color_by == "source" and "source" not in table.columns:
        raise ValueError("color_by='source' needs a `source` column in the metrics table")

    logger.info(
        f"Loaded {len(table):,} groups; biology panels: "
        f"{[column for column, _, _ in biology_panels(table)] or 'none'}"
    )
    if config.color_by == "source":
        logger.info(f"Comparison mode: {len(_source_series(table))} series {_source_series(table)}")
    plot_coherence(table, config.labels, config.color_by)
    save_dual(config.output.with_suffix(""))
    logger.success(f"Wrote {config.output}")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Generate coherence visualization")
    parser.add_argument("--input", dest="input_metrics", type=Path, required=True,
                        help="Input coherence metrics table (per-source Parquet, or the cross-source combined Parquet)")
    parser.add_argument("--color-by", dest="color_by", choices=("none", "source"), default="none",
                        help="'none' = the single-series figure (z as colour); 'source' = the cross-source comparison, one colour per source")
    parser.add_argument("--dedup-series", dest="dedup_series", type=Path, default=None,
                        help="Optional de-duplicated representatives TSV, appended as a `dedup` series (comparison mode)")
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
            color_by=args.color_by,
            dedup_series=args.dedup_series,
        )
        run(config)
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
