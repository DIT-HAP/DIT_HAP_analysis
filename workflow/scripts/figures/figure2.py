#!/usr/bin/env python3

"""
Figure 2 — Depletion kinetics: how one fitted curve defines DR, DL and A
========================================================================

Panel A is the model-definition panel that every downstream DR/DL feature is
read from. It has two halves joined by an arrow:

  left   the analysis flow — time-series LFC data, the Gompertz equation with
         the three parameter names colour-coded, and the extracted parameters;
  right  one example fit drawn from the release: the measured points and the
         fitted curve, with the three parameters written onto the geometry that
         defines them — A as the plateau line, DR as the slope of the tangent
         at the inflection, DL as that tangent's x-intercept.

Panel B sits directly below A at the same width: four example fits (tna1,
pst2, rpl3801, nrs1) in a 2 x 2 grid, each mini axes repeating panel A's plot
footprint. The grid shares one x scale and one y scale so the four curves are
directly comparable; only the left column carries y tick labels and only the
bottom row carries x tick labels. Each mini axes repeats the three guides in a
lighter weight; their text labels stay in A, and the four genes' fitted values
and R2 share one table under the grid.

Panel C sits to the right of A on the same row: the distributions of DR and DL
across every successful fit in the release, one histogram each, sharing one
count scale so the two shapes are directly comparable. It answers "where in this
space do real genes fall?" — the spread A's single example cannot show.

The genes are fixed (``PANEL_A_GENE`` and ``PANEL_B_GENES``): they are the
instances the figure tells its story with, not per-run choices. Each curve is
evaluated from its released amplitude, rate and lag (A, DR, DL); the five
``*_fitted`` columns are used only to confirm this script's model reproduces
the fit of record. Every drawing decision lives in this file; the script
imports no stage module from ``workflow/src`` and calls no upstream renderer.

The parameter colours are one system across the figure: each parameter has a
palette colour, its name in the equation sits on a tint of that colour, and the
matching guide in the plot (line, dots, label) uses it. Data and fit use
neutral colours so the three parameter colours never compete with the data.

The page is an A4 manuscript page (210 x 297 mm). ``cns.multipanel`` sizes only
the canvas width and lets the rows decide the height, so the flow below draws a
content-sized canvas and the assembly mounts it 40 pt in from the A4 page's
top-left corner, leaving the rest of the page blank. The panel letter is placed
by multipanel at the left edge of the panel's column.

Input
-----
- <HD_DIT_HAP release>/gene_level/gene_level_fitting_statistics.tsv
  Columns: Name, time_points, A, DR, DL, R2, t_inflection, y_inflection,
  YES0..YES4, YES0_fitted..YES4_fitted (panel C additionally reads Status).

Output
------
- results/7a_figures/figure2.pdf                content mounted on an A4 page
- results/7a_figures/figure2.review.png         360 dpi review raster
- results/7a_figures/_work/figure2_canvas.pdf   canvas before the A4 mount (debug)

Usage
-----
    python workflow/scripts/figures/figure2.py \\
        --gene-fitting <release>/gene_level/gene_level_fitting_statistics.tsv \\
        --output-pdf results/7a_figures/figure2.pdf \\
        --output-review results/7a_figures/figure2.review.png \\
        --canvas results/7a_figures/_work/figure2_canvas.pdf

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-13
Version:  2.0.0
"""

# =============================================================================
# 1. IMPORTS (Standard Library -> Third-Party)
# =============================================================================
import argparse
from dataclasses import dataclass
from math import ceil
from pathlib import Path
import subprocess
import sys
from typing import cast

import cnsplots as cns
from loguru import logger
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_rgba
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.text import Text
import numpy as np
import pandas as pd
from pypdf import PdfReader, PdfWriter, Transformation

# =============================================================================
# 2. CONSTANTS (house style, layout, draw vocabulary)
# =============================================================================
# --- House style (self-contained copy; see docs/figures.md) -----------------
HOUSE_PALETTE = "Cell"
SEQUENTIAL_CMAP = "viridis"
FONT_SANS_SERIF = (
    "Arial",
    "Helvetica",
    "Nimbus Sans",
    "Liberation Sans",
    "DejaVu Sans",
)
HOUSE_FIGURE_DPI = 360
FURNITURE_COLOR = cns.GRAY  # box edges, arrows, guides that carry no parameter

# --- Canvas and panel geometry, in points -----------------------------------
# The manuscript page is A4 (210 x 297 mm). multipanel sizes the figure by
# width only, so the content canvas is A4 wide minus the page margins and the
# assembly mounts it PAGE_MARGIN_PT in from the page's top-left corner; the
# rest of the page stays blank.
MM_PER_INCH = 25.4
A4_WIDTH_PT = 210.0 / MM_PER_INCH * 72.0  # 595.28
A4_HEIGHT_PT = 297.0 / MM_PER_INCH * 72.0  # 841.89
PAGE_MARGIN_PT = 40.0  # white border at the left, top and right
CANVAS_WIDTH_PT = A4_WIDTH_PT - 2 * PAGE_MARGIN_PT  # 515.28

# Panel A is one container axes holding the whole composition. Its width is
# the flow column (80) + the arrow (20) + the plot's y-axis decoration (~24)
# + the plot (100) + a 6 pt right margin.
PANEL_WIDTH_PT = 236.0
PANEL_HEIGHT_PT = 100.0

# Left flow column: three boxes stacked, joined by short vertical arrows.
FLOW_X_PT = 2.0
FLOW_WIDTH_PT = 78.0
FLOW_CENTER_X_PT = FLOW_X_PT + FLOW_WIDTH_PT / 2
DATA_BOX_Y_PT = 80.0
DATA_BOX_HEIGHT_PT = 17.0
FORMULA_BOX_Y_PT = 40.0
FORMULA_BOX_HEIGHT_PT = 33.5
PARAMETER_BOX_Y_PT = 14.0
PARAMETER_BOX_HEIGHT_PT = 19.5
BOX_FACE_COLOR = "#f7f7f7"
BOX_EDGE_PT = 0.6
ARROW_HEAD_SCALE = 6.0
ARROW_LINE_WIDTH_PT = 0.8
# Horizontal arrow from the flow to the plot; 20 pt is the gap the composition
# was laid out around.
FLOW_ARROW_Y_PT = 56.0
FLOW_ARROW_END_X_PT = 101.0
FLOW_ARROW_HEAD_SCALE = 8.0

# Right plot: a 100 x 76 pt axes inset into the container. The container owns
# the plot's y-axis decoration and x label, so the panel box is exactly
# PANEL_WIDTH_PT x PANEL_HEIGHT_PT including them.
PLOT_X_PT = 126.0
PLOT_Y_PT = 22.0
PLOT_WIDTH_PT = 100.0
PLOT_HEIGHT_PT = 76.0

# Panel B sits below A at the same width: four example fits in a 2 x 2 grid of
# mini axes, each repeating panel A's plot footprint (100 x 76 pt). Only the
# left column carries y tick labels and only the bottom row x tick labels, so
# the two 100 pt columns fit inside the same 236 pt width.
PANEL_B_WIDTH_PT = PANEL_WIDTH_PT
B_GRID_COLUMNS = 2
B_MINI_WIDTH_PT = PLOT_WIDTH_PT
B_MINI_HEIGHT_PT = PLOT_HEIGHT_PT
B_LEFT_DECORATION_PT = 25.0
B_COLUMN_GAP_PT = 6.0
B_ROW_GAP_PT = 14.0  # clears the lower row's gene title
B_TOP_TITLE_BAND_PT = 12.0
B_BOTTOM_DECORATION_PT = 22.0
B_Y_LABEL_CENTER_X_PT = 12.0
B_X_LABEL_BOTTOM_PT = 3.0
B_GUIDE_LINE_WIDTH_PT = 0.8
B_GUIDE_DOT_SIZE = 8.0
B_Y_PADDING_BOTTOM = 0.08
B_Y_PADDING_TOP = 0.12
# The four genes' fitted values live in one table under the grid instead of
# four per-axes blocks; the band holds the table plus the gap above it.
B_TABLE_BAND_PT = 58.0
B_TABLE_BOTTOM_PT = 2.0
B_TABLE_HEADER_HEIGHT_PT = 10.0
B_TABLE_ROW_HEIGHT_PT = 9.0
B_TABLE_FONT_PT = 6.5
B_TABLE_GENE_X_PT = 42.0
B_TABLE_AMPLITUDE_RIGHT_PT = 111.0
B_TABLE_RATE_RIGHT_PT = 144.0
B_TABLE_LAG_RIGHT_PT = 177.0
B_TABLE_R2_RIGHT_PT = 214.0
B_TABLE_RULE_Y_PT = 38.0
PANEL_B_HEIGHT_PT = (
    B_TABLE_BAND_PT
    + B_BOTTOM_DECORATION_PT
    + 2 * B_MINI_HEIGHT_PT
    + B_ROW_GAP_PT
    + B_TOP_TITLE_BAND_PT
)
# Panel C sits to the right of A on the same row: two histograms, DR and DL,
# each at panel A's plot footprint. The two panels have independent count
# axes because DL's zero spike (3,113 genes) is far taller than DR's peak.
PANEL_C_WIDTH_PT = 252.0
PANEL_C_HEIGHT_PT = 100.0
C_LEFT_DECORATION_PT = 34.0  # count labels are wider than LFC ticks + ylabel
C_COLUMN_GAP_PT = 10.0
C_PLOT_Y_PT = 22.0
C_Y_LABEL_CENTER_X_PT = 10.0
C_Y_HEADROOM = 1.08

ROW_GAP_PT = 10.0  # between panel A and panel B

# --- Panel vocabulary -------------------------------------------------------
PANEL_A_GENE = "mtf1"
PANEL_B_GENES = ("tna1", "pst2", "rpl3801", "nrs1")

# --- Panel C vocabulary -----------------------------------------------------
DISTRIBUTION_COLUMNS = ("Status", "DR", "DL")
SUCCESS_STATUS = "Success"
DR_HIST_RANGE = (-2.2, 0.3)
DR_BIN_WIDTH = 0.1
DL_HIST_RANGE = (0.0, 10.0)
DL_BIN_WIDTH = 0.5
DR_X_LABEL = "Depletion rate (DR)"
DL_X_LABEL = "Depletion lag (DL)"
COUNT_Y_LABEL = "Number of genes"
NAME_COLUMN = "Name"
TIME_POINTS_COLUMN = "time_points"
AMPLITUDE_COLUMN = "A"
RATE_COLUMN = "DR"
LAG_COLUMN = "DL"
R2_COLUMN = "R2"
T_INFLECTION_COLUMN = "t_inflection"
Y_INFLECTION_COLUMN = "y_inflection"
VALUE_COLUMNS = ("YES0", "YES1", "YES2", "YES3", "YES4")
FITTED_COLUMNS = tuple(f"{column}_fitted" for column in VALUE_COLUMNS)
# The release keeps only a few decimals for A/DR/DL and the fitted values, so
# evaluating the model from the stored parameters can sit a few units in the
# last place away from the stored fit; 5e-3 covers the rounding of every
# parameter without accepting a different model.
FITTED_TOLERANCE = 5e-3

# --- Parameter colours ------------------------------------------------------
# One colour per parameter, used by the equation chip, the guide and the label.
# Cell slots 2/0/1 are the greyscale-separable ochre-red-teal triad.
AMPLITUDE_COLOR = cns.get_hexcolors_from_apalette([2], HOUSE_PALETTE)[0]
RATE_COLOR = cns.get_hexcolors_from_apalette([0], HOUSE_PALETTE)[0]
LAG_COLOR = cns.get_hexcolors_from_apalette([1], HOUSE_PALETTE)[0]
# Data and fit take neutral slots so no parameter colour is a data colour.
FIT_COLOR = cns.get_hexcolors_from_apalette([3], HOUSE_PALETTE)[0]
POINT_COLOR = cns.get_hexcolors_from_apalette([9], HOUSE_PALETTE)[0]
CHIP_ALPHA = 0.25

# --- Panel A drawing sizes --------------------------------------------------
FLOW_TEXT_FONT_PT = 6.5
EQUATION_FONT_PT = 6.5
PLOT_TICK_FONT_PT = 6.5
PLOT_LABEL_FONT_PT = 7.0
EQUATION_LINE_GAP_PT = 0.7
CHIP_PAD_PT = 0.8
CHIP_ROUNDING_PT = 1.0
CURVE_SAMPLES = 300
OBSERVED_MARKER_SIZE = 12.0
OBSERVED_MARKER_EDGE_PT = 0.4
GUIDE_LINE_WIDTH_PT = 1.0
GUIDE_DOT_SIZE = 16.0
ZERO_LINE_WIDTH_PT = 0.5

# --- Panel A text -----------------------------------------------------------
DATA_BOX_TEXT = "Time-series LFC data"
PARAMETER_BOX_TEXT = "Extracted parameters"
PARAMETER_NAMES = (
    (r"$DR$", RATE_COLOR),
    (r"$DL$", LAG_COLOR),
    (r"$A$", AMPLITUDE_COLOR),
)
PARAMETER_NAME_FONT_PT = 7.0
PARAMETER_NAME_GAP_PT = 9.0
X_LABEL = "Generations"
Y_LABEL = r"LFC ($\log_2$)"
PLATEAU_LABEL = r"$A$ (plateau)"
RATE_LABEL = r"$DR$ (slope)"
LAG_LABEL = r"$DL$ (intercept)"
# One decimal is the display precision the release's generation grid needs
# (2.352 -> 2.4); a trailing ".0" is dropped so the first tick reads 0.
TICK_DECIMALS = 1
Y_TICK_STEP = 2.0
# The data fill the middle band; the top is left for the tangent/labels and the
# bottom for the plateau label and the delay line's lower dot.
Y_PADDING_BOTTOM = 0.28
Y_PADDING_TOP = 0.31
X_PADDING = 0.03
# Panel B needs more room on the left than A: three of its four genes have
# DL = 0, and the delay guide must stay visibly inside the axes instead of
# hugging the spine.
B_X_PADDING_LEFT = 0.09
B_X_PADDING_RIGHT = 0.03
DELAY_BOTTOM_MARGIN = 0.15  # keeps the delay line's lower dot inside the axes

# --- Input contract ---------------------------------------------------------
GENE_FIT_COLUMNS = (
    NAME_COLUMN,
    TIME_POINTS_COLUMN,
    AMPLITUDE_COLUMN,
    RATE_COLUMN,
    LAG_COLUMN,
    R2_COLUMN,
    T_INFLECTION_COLUMN,
    Y_INFLECTION_COLUMN,
    *VALUE_COLUMNS,
    *FITTED_COLUMNS,
)


# =============================================================================
# 3. CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class Figure2Config:
    """Resolved input/output paths for Figure 2."""

    gene_fitting: Path
    output_pdf: Path
    output_review: Path
    canvas_pdf: Path

    def validate(self) -> None:
        """Raise if a required input is missing, then create the output dirs."""
        if not self.gene_fitting.exists():
            raise ValueError(f"Required input not found: {self.gene_fitting}")
        for output in [self.output_pdf, self.output_review, self.canvas_pdf]:
            output.parent.mkdir(parents=True, exist_ok=True)


@dataclass(kw_only=True, slots=True, frozen=True)
class GeneFit:
    """One gene's measured time course and its released Gompertz parameters."""

    name: str
    generations: np.ndarray
    observed: np.ndarray
    amplitude: float
    rate: float
    lag: float
    r2: float
    t_inflection: float
    y_inflection: float


@dataclass(frozen=True)
class FormulaPiece:
    """One run of text in the equation, optionally sitting on a colour chip."""

    text: str
    color: str = "black"
    chip: tuple[float, float, float, float] | None = None


# =============================================================================
# 4. DATA PREPARATION (pure; schemas are validated at the I/O boundary)
# =============================================================================
def gompertz(generations: np.ndarray, amplitude: float, rate: float, lag: float) -> np.ndarray:
    """Evaluate the release's Gompertz model, clipped so exp cannot overflow."""
    if amplitude == 0:
        return np.zeros_like(generations)
    alpha = (rate * np.e) / amplitude
    exponent = alpha * (lag - generations) + 1
    return amplitude * np.exp(-np.exp(np.clip(exponent, -700, 700)))


def build_gene_fit(fitting: pd.DataFrame, gene_name: str) -> GeneFit:
    """Select one gene and tie this script's model to its released fit."""
    is_example = fitting[NAME_COLUMN] == gene_name
    matches = fitting.loc[is_example]
    if matches.empty:
        raise ValueError(f"{gene_name} is not in the gene-level fitting table")
    if len(matches) > 1:
        raise ValueError(f"{gene_name} appears {len(matches)} times in the fitting table")

    row = matches.iloc[0]
    generations = np.array(
        [float(point) for point in row[TIME_POINTS_COLUMN].split(",")], dtype=float
    )
    if len(generations) != len(VALUE_COLUMNS):
        raise ValueError(
            f"{gene_name} has {len(generations)} time points for {len(VALUE_COLUMNS)} values"
        )

    amplitude = float(row[AMPLITUDE_COLUMN])
    rate = float(row[RATE_COLUMN])
    lag = float(row[LAG_COLUMN])
    released_fitted = row[list(FITTED_COLUMNS)].to_numpy(dtype=float)
    model_fitted = gompertz(generations, amplitude, rate, lag)
    deviation = float(np.max(np.abs(model_fitted - released_fitted)))
    if deviation > FITTED_TOLERANCE:
        raise ValueError(
            f"{gene_name}: model from A/DR/DL deviates {deviation:.4f} from the released fit"
        )

    return GeneFit(
        name=str(row[NAME_COLUMN]),
        generations=generations,
        observed=row[list(VALUE_COLUMNS)].to_numpy(dtype=float),
        amplitude=amplitude,
        rate=rate,
        lag=lag,
        r2=float(row[R2_COLUMN]),
        t_inflection=float(row[T_INFLECTION_COLUMN]),
        y_inflection=float(row[Y_INFLECTION_COLUMN]),
    )


def build_parameter_distributions(table: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Return DR and DL across every successful gene-level fit."""
    is_success = table["Status"] == SUCCESS_STATUS
    selected = table.loc[is_success]
    return selected["DR"].to_numpy(dtype=float), selected["DL"].to_numpy(dtype=float)


def histogram_edges(low: float, high: float, width: float) -> np.ndarray:
    """Return bin edges covering [low, high], high included in the last bin."""
    return np.arange(low, high + width / 2, width)


def display_ticks(generations: np.ndarray) -> list[str]:
    """Format the generation grid at the panel's display precision."""
    labels = []
    for generation in generations:
        label = f"{generation:.{TICK_DECIMALS}f}"
        labels.append(label.removesuffix(".0"))
    return labels


def y_data_bounds(fit: GeneFit, curve: np.ndarray) -> tuple[float, float]:
    """Return the lowest and highest y among the points, the curve and the plateau."""
    # The plateau guide is drawn at A even when the sampled curve has not
    # reached it by the last generation, so A must be inside the axes.
    low = min(float(fit.observed.min()), float(curve.min()), fit.amplitude)
    high = max(float(fit.observed.max()), float(curve.max()))
    return low, high


def padded_y_limits(low: float, high: float) -> tuple[float, float]:
    """Pad the data bounds, leaving the wider band on top for the tangent."""
    span = high - low
    return low - Y_PADDING_BOTTOM * span, high + Y_PADDING_TOP * span


def y_tick_values(low: float, high: float) -> np.ndarray:
    """Tick the 2-unit grid that falls inside the data bounds."""
    first = np.ceil(low / Y_TICK_STEP) * Y_TICK_STEP
    last = np.floor(high / Y_TICK_STEP) * Y_TICK_STEP
    return np.arange(first, last + Y_TICK_STEP / 2, Y_TICK_STEP)


def equation_lines() -> tuple[tuple[FormulaPiece, ...], tuple[FormulaPiece, ...]]:
    """Return the Gompertz equation split over two lines, the parameters chipped."""
    line_one = (
        FormulaPiece(r"$\mathrm{LFC}(t)=$"),
        FormulaPiece(r"$A$", AMPLITUDE_COLOR, to_rgba(AMPLITUDE_COLOR, CHIP_ALPHA)),
        FormulaPiece(r"$\cdot\exp(-\exp[$"),
    )
    line_two = (
        FormulaPiece(r"$($"),
        FormulaPiece(r"$DR$", RATE_COLOR, to_rgba(RATE_COLOR, CHIP_ALPHA)),
        FormulaPiece(r"$\cdot e/$"),
        FormulaPiece(r"$A$", AMPLITUDE_COLOR, to_rgba(AMPLITUDE_COLOR, CHIP_ALPHA)),
        FormulaPiece(r"$)($"),
        FormulaPiece(r"$DL$", LAG_COLOR, to_rgba(LAG_COLOR, CHIP_ALPHA)),
        FormulaPiece(r"$-t)+1])$"),
    )
    return line_one, line_two


# =============================================================================
# 5. RENDERING
# =============================================================================
def apply_house_style() -> None:
    """Configure cnsplots with the project house style (Cell palette, Arial first)."""
    cns.settings.font_sans_serif = FONT_SANS_SERIF
    cns.settings.panel_label_fontname = None
    cns.settings.pdf_fonttype = 42  # TrueType embedding, editable in Illustrator
    cns.settings.svg_fonttype = "none"
    cns.settings.figure_dpi = HOUSE_FIGURE_DPI
    cns.settings.savefig_dpi = HOUSE_FIGURE_DPI
    cns.settings.savefig_transparent = False
    cns.settings.palette_qual = HOUSE_PALETTE
    cns.settings.palette_seq = SEQUENTIAL_CMAP
    cns.setup_matplotlib(color_cycle=HOUSE_PALETTE, color_map=SEQUENTIAL_CMAP)


def draw_figure(
    panel_a_fit: GeneFit,
    panel_b_fits: tuple[GeneFit, ...],
    distributions: tuple[np.ndarray, np.ndarray],
) -> cns.multipanel:
    """Draw panels A and B and return the multipanel that owns the page.

    Returns the multipanel itself, not just its figure: cnsplots keeps the
    measured-layout draw callback as a weak reference, so the caller must hold
    the object until the file is written or the panels stay unmeasured.
    """
    mp = cns.multipanel(max_width=CANVAS_WIDTH_PT, unit="pt")
    container = mp.panel("A", width=PANEL_WIDTH_PT, height=PANEL_HEIGHT_PT, margin_bottom=0)
    blank_container(container)
    render_panel(container, panel_a_fit)

    container_c = mp.panel(
        "C",
        width=PANEL_C_WIDTH_PT,
        height=PANEL_C_HEIGHT_PT,
        margin_right=0,
        margin_bottom=0,
    )
    blank_container(container_c)
    render_panel_c(container_c, distributions)

    mp.newline()
    container_b = mp.panel(
        "B",
        width=PANEL_B_WIDTH_PT,
        height=PANEL_B_HEIGHT_PT,
        margin_top=ROW_GAP_PT,
        margin_bottom=0,
    )
    blank_container(container_b)
    render_panel_b(container_b, panel_b_fits)
    return mp


def blank_container(container: Axes) -> None:
    """Clear the ticks, turn the axis off and give the container a pt coordinate system."""
    # Clear the ticks before turning the axis off: multipanel measures y-axis
    # decoration even on an axis-off axes, and placeholder ticks would widen
    # the panel's left reserve.
    container.set_xticks([])
    container.set_yticks([])
    container.set_axis_off()


def render_panel(container: Axes, fit: GeneFit) -> None:
    """Draw the flow column, the horizontal arrow and the anatomy plot."""
    # The container's data coordinates are points, so every box and arrow below
    # is placed in the panel's own pt geometry.
    container.set_xlim(0, PANEL_WIDTH_PT)
    container.set_ylim(0, PANEL_HEIGHT_PT)
    draw_flow(container, fit)
    draw_horizontal_arrow(container)

    plot = container.inset_axes(
        (
            PLOT_X_PT / PANEL_WIDTH_PT,
            PLOT_Y_PT / PANEL_HEIGHT_PT,
            PLOT_WIDTH_PT / PANEL_WIDTH_PT,
            PLOT_HEIGHT_PT / PANEL_HEIGHT_PT,
        )
    )
    render_anatomy_plot(plot, fit)


def render_panel_c(container: Axes, distributions: tuple[np.ndarray, np.ndarray]) -> None:
    """Draw the DR and DL histograms side by side at panel A's plot footprint."""
    container.set_xlim(0, PANEL_C_WIDTH_PT)
    container.set_ylim(0, PANEL_C_HEIGHT_PT)
    dr_values, dl_values = distributions
    panels = (
        (dr_values, DR_HIST_RANGE, DR_BIN_WIDTH, RATE_COLOR, DR_X_LABEL),
        (dl_values, DL_HIST_RANGE, DL_BIN_WIDTH, LAG_COLOR, DL_X_LABEL),
    )
    edges = [
        histogram_edges(hist_range[0], hist_range[1], bin_width)
        for _, hist_range, bin_width, _, _ in panels
    ]
    counts = [
        np.histogram(values, bins=bins)[0]
        for values, bins in zip((dr_values, dl_values), edges, strict=True)
    ]
    # One count scale for both panels: the two distributions describe the same
    # 4,513 genes, and the shared scale makes DL's zero spike comparable to DR.
    y_max = max(int(count.max()) for count in counts) * C_Y_HEADROOM
    for column, ((values, _, _, color, xlabel), bins) in enumerate(zip(panels, edges, strict=True)):
        x_pt = C_LEFT_DECORATION_PT + column * (PLOT_WIDTH_PT + C_COLUMN_GAP_PT)
        ax = container.inset_axes(
            (
                x_pt / PANEL_C_WIDTH_PT,
                C_PLOT_Y_PT / PANEL_C_HEIGHT_PT,
                PLOT_WIDTH_PT / PANEL_C_WIDTH_PT,
                PLOT_HEIGHT_PT / PANEL_C_HEIGHT_PT,
            )
        )
        render_distribution_axes(ax, values, bins, color, xlabel, y_max, show_y_labels=column == 0)
    container.text(
        C_Y_LABEL_CENTER_X_PT,
        C_PLOT_Y_PT + PLOT_HEIGHT_PT / 2,
        COUNT_Y_LABEL,
        rotation=90,
        ha="center",
        va="center",
        fontsize=PLOT_LABEL_FONT_PT,
    )


def render_distribution_axes(
    ax: Axes,
    values: np.ndarray,
    edges: np.ndarray,
    color: str,
    xlabel: str,
    y_max: float,
    show_y_labels: bool,
) -> None:
    """Draw one parameter's count histogram under the shared count scale."""
    cns.histplot(x=values, bins=edges, color=color, ax=ax, legend=False)
    ax.set_ylim(0, y_max)
    ax.set_ylabel("")
    ax.set_xlabel(xlabel, fontsize=PLOT_LABEL_FONT_PT)
    ax.tick_params(axis="both", labelsize=PLOT_TICK_FONT_PT, labelleft=show_y_labels)


def render_panel_b(container: Axes, fits: tuple[GeneFit, ...]) -> None:
    """Draw the 2 x 2 grid of example fits at panel A's plot footprint."""
    container.set_xlim(0, PANEL_B_WIDTH_PT)
    container.set_ylim(0, PANEL_B_HEIGHT_PT)
    n_rows = ceil(len(fits) / B_GRID_COLUMNS)
    data_bounds, y_limits = panel_b_y_limits(fits)
    x_limits = panel_b_x_limits(fits[0])
    for index, fit in enumerate(fits):
        row, column = divmod(index, B_GRID_COLUMNS)
        x_pt, y_pt = panel_b_position(index, n_rows)
        ax = container.inset_axes(
            (
                x_pt / PANEL_B_WIDTH_PT,
                y_pt / PANEL_B_HEIGHT_PT,
                B_MINI_WIDTH_PT / PANEL_B_WIDTH_PT,
                B_MINI_HEIGHT_PT / PANEL_B_HEIGHT_PT,
            )
        )
        render_example_axes(
            ax,
            fit,
            x_limits,
            y_limits,
            data_bounds,
            show_y_labels=column == 0,
            show_x_labels=row == n_rows - 1,
        )
    add_panel_b_labels(container, n_rows)
    add_panel_b_table(container, fits)


def panel_b_position(index: int, n_rows: int) -> tuple[float, float]:
    """Return one mini axes' lower-left corner in panel B points, top row first."""
    row, column = divmod(index, B_GRID_COLUMNS)
    x_pt = B_LEFT_DECORATION_PT + column * (B_MINI_WIDTH_PT + B_COLUMN_GAP_PT)
    y_pt = (
        B_TABLE_BAND_PT
        + B_BOTTOM_DECORATION_PT
        + (n_rows - 1 - row) * (B_MINI_HEIGHT_PT + B_ROW_GAP_PT)
    )
    return x_pt, y_pt


def panel_b_y_limits(
    fits: tuple[GeneFit, ...],
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Return the shared data bounds and the padded y limits across the example fits."""
    bounds = [y_data_bounds(fit, fit_curve(fit)[1]) for fit in fits]
    low = min(bound[0] for bound in bounds)
    high = max(bound[1] for bound in bounds)
    span = high - low
    return (low, high), (low - B_Y_PADDING_BOTTOM * span, high + B_Y_PADDING_TOP * span)


def render_example_axes(
    ax: Axes,
    fit: GeneFit,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
    data_bounds: tuple[float, float],
    show_y_labels: bool,
    show_x_labels: bool,
) -> None:
    """Draw one example gene under the grid's shared axes conventions."""
    curve_x, curve_y = fit_curve(fit)
    draw_fit_series(ax, fit, curve_x, curve_y)
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    draw_plateau_guide(ax, fit, x_limits[0], x_limits[1], with_label=False)
    draw_rate_guide(ax, fit, x_limits[0], x_limits[1], y_limits[1], y_limits[0], with_label=False)
    draw_delay_guide(ax, fit, y_limits[0], with_label=False)
    ax.set_xticks(fit.generations, display_ticks(fit.generations))
    ax.set_yticks(y_tick_values(*data_bounds))
    ax.tick_params(
        axis="both",
        labelsize=PLOT_TICK_FONT_PT,
        labelleft=show_y_labels,
        labelbottom=show_x_labels,
    )
    ax.set_title(fit.name, fontstyle="italic", fontweight="bold")


def add_panel_b_labels(container: Axes, n_rows: int) -> None:
    """Centre one x label and one rotated y label on the whole grid."""
    grid_height_pt = n_rows * B_MINI_HEIGHT_PT + (n_rows - 1) * B_ROW_GAP_PT
    grid_width_pt = B_GRID_COLUMNS * B_MINI_WIDTH_PT + (B_GRID_COLUMNS - 1) * B_COLUMN_GAP_PT
    container.text(
        B_Y_LABEL_CENTER_X_PT,
        B_TABLE_BAND_PT + B_BOTTOM_DECORATION_PT + grid_height_pt / 2,
        Y_LABEL,
        rotation=90,
        ha="center",
        va="center",
        fontsize=PLOT_LABEL_FONT_PT,
    )
    container.text(
        B_LEFT_DECORATION_PT + grid_width_pt / 2,
        B_TABLE_BAND_PT + B_X_LABEL_BOTTOM_PT,
        X_LABEL,
        ha="center",
        va="bottom",
        fontsize=PLOT_LABEL_FONT_PT,
    )


def add_panel_b_table(container: Axes, fits: tuple[GeneFit, ...]) -> None:
    """Write the four genes' fitted values and R2 in one table under the grid."""
    table_height = B_TABLE_HEADER_HEIGHT_PT + len(fits) * B_TABLE_ROW_HEIGHT_PT
    header_y = B_TABLE_BOTTOM_PT + table_height - B_TABLE_HEADER_HEIGHT_PT / 2
    container.text(
        B_TABLE_GENE_X_PT,
        header_y,
        "Gene",
        ha="left",
        va="center",
        fontsize=B_TABLE_FONT_PT,
        color=FURNITURE_COLOR,
    )
    headers = (
        (r"$A$", B_TABLE_AMPLITUDE_RIGHT_PT, AMPLITUDE_COLOR),
        (r"$DR$", B_TABLE_RATE_RIGHT_PT, RATE_COLOR),
        (r"$DL$", B_TABLE_LAG_RIGHT_PT, LAG_COLOR),
        (r"$R^2$", B_TABLE_R2_RIGHT_PT, FURNITURE_COLOR),
    )
    for text, x_pt, color in headers:
        container.text(
            x_pt, header_y, text, ha="right", va="center", fontsize=B_TABLE_FONT_PT, color=color
        )
    container.plot(
        [B_TABLE_GENE_X_PT, B_TABLE_R2_RIGHT_PT],
        [B_TABLE_RULE_Y_PT, B_TABLE_RULE_Y_PT],
        color=FURNITURE_COLOR,
        linewidth=0.4,
        zorder=1,
    )

    for index, fit in enumerate(fits):
        row_y = (
            header_y
            - B_TABLE_HEADER_HEIGHT_PT / 2
            - B_TABLE_ROW_HEIGHT_PT / 2
            - index * B_TABLE_ROW_HEIGHT_PT
        )
        lag = 0.0 if fit.lag == 0 else fit.lag
        container.text(
            B_TABLE_GENE_X_PT,
            row_y,
            fit.name,
            ha="left",
            va="center",
            fontsize=B_TABLE_FONT_PT,
            fontstyle="italic",
        )
        values = (
            (f"{fit.amplitude:.2f}", B_TABLE_AMPLITUDE_RIGHT_PT, AMPLITUDE_COLOR),
            (f"{fit.rate:.2f}", B_TABLE_RATE_RIGHT_PT, RATE_COLOR),
            (f"{lag:.2f}", B_TABLE_LAG_RIGHT_PT, LAG_COLOR),
            (f"{fit.r2:.3f}", B_TABLE_R2_RIGHT_PT, FURNITURE_COLOR),
        )
        for text, x_pt, color in values:
            container.text(
                x_pt, row_y, text, ha="right", va="center", fontsize=B_TABLE_FONT_PT, color=color
            )


def draw_flow(container: Axes, fit: GeneFit) -> None:
    """Draw the three stacked flow boxes, their arrows, and the fitted equation."""
    draw_box(container, DATA_BOX_Y_PT, DATA_BOX_HEIGHT_PT)
    container.text(
        FLOW_CENTER_X_PT,
        DATA_BOX_Y_PT + DATA_BOX_HEIGHT_PT / 2,
        DATA_BOX_TEXT,
        ha="center",
        va="center",
        fontsize=FLOW_TEXT_FONT_PT,
    )
    draw_vertical_arrow(container, FORMULA_BOX_Y_PT + FORMULA_BOX_HEIGHT_PT, DATA_BOX_Y_PT)

    draw_box(container, FORMULA_BOX_Y_PT, FORMULA_BOX_HEIGHT_PT)
    draw_equation(container, fit)

    draw_vertical_arrow(container, PARAMETER_BOX_Y_PT + PARAMETER_BOX_HEIGHT_PT, FORMULA_BOX_Y_PT)
    draw_box(container, PARAMETER_BOX_Y_PT, PARAMETER_BOX_HEIGHT_PT)
    container.text(
        FLOW_CENTER_X_PT,
        PARAMETER_BOX_Y_PT + PARAMETER_BOX_HEIGHT_PT - 4.5,
        PARAMETER_BOX_TEXT,
        ha="center",
        va="center",
        fontsize=FLOW_TEXT_FONT_PT,
    )
    draw_parameter_names(container)


def draw_box(container: Axes, y_pt: float, height_pt: float) -> None:
    """Draw one rounded flow box in the left column."""
    container.add_patch(
        FancyBboxPatch(
            (FLOW_X_PT, y_pt),
            FLOW_WIDTH_PT,
            height_pt,
            boxstyle="round,pad=0,rounding_size=2.0",
            linewidth=BOX_EDGE_PT,
            edgecolor=FURNITURE_COLOR,
            facecolor=BOX_FACE_COLOR,
            zorder=1,
        )
    )


def draw_vertical_arrow(container: Axes, tail_pt: float, head_pt: float) -> None:
    """Draw the arrow between two stacked flow boxes."""
    container.add_patch(
        FancyArrowPatch(
            (FLOW_CENTER_X_PT, tail_pt),
            (FLOW_CENTER_X_PT, head_pt),
            arrowstyle="-|>",
            mutation_scale=ARROW_HEAD_SCALE,
            color=FURNITURE_COLOR,
            linewidth=ARROW_LINE_WIDTH_PT,
            shrinkA=0,
            shrinkB=0,
            zorder=2,
        )
    )


def draw_horizontal_arrow(container: Axes) -> None:
    """Draw the arrow from the flow column to the anatomy plot."""
    container.add_patch(
        FancyArrowPatch(
            (FLOW_X_PT + FLOW_WIDTH_PT, FLOW_ARROW_Y_PT),
            (FLOW_ARROW_END_X_PT, FLOW_ARROW_Y_PT),
            arrowstyle="-|>",
            mutation_scale=FLOW_ARROW_HEAD_SCALE,
            color=FURNITURE_COLOR,
            linewidth=ARROW_LINE_WIDTH_PT,
            shrinkA=0,
            shrinkB=0,
            zorder=2,
        )
    )


def text_extent_pt(container: Axes, text: str, fontsize: float) -> tuple[float, float]:
    """Measure a text's rendered width and height in points."""
    renderer = cast(FigureCanvasAgg, container.figure.canvas).get_renderer()
    artist = container.text(0, 0, text, fontsize=fontsize)
    bbox = artist.get_window_extent(renderer)
    artist.remove()
    return bbox.width * 72 / container.figure.dpi, bbox.height * 72 / container.figure.dpi


def draw_equation(container: Axes, fit: GeneFit) -> None:
    """Draw the two equation lines, centred in the formula box, with colour chips."""
    line_one, line_two = equation_lines()
    line_one_y = FORMULA_BOX_Y_PT + 25.3
    line_two_y = FORMULA_BOX_Y_PT + 16.3
    draw_equation_line(container, line_one, line_one_y)
    draw_equation_line(container, line_two, line_two_y)
    container.text(
        FLOW_CENTER_X_PT,
        FORMULA_BOX_Y_PT + 8.5,
        rf"$R^2={fit.r2:.3f}$",
        ha="center",
        va="center",
        fontsize=FLOW_TEXT_FONT_PT,
        color=FURNITURE_COLOR,
    )


def draw_equation_line(container: Axes, pieces: tuple[FormulaPiece, ...], y_pt: float) -> None:
    """Lay one equation line out left-to-right, centred, drawing each chip behind its text."""
    widths = [text_extent_pt(container, piece.text, EQUATION_FONT_PT)[0] for piece in pieces]
    total = sum(widths) + EQUATION_LINE_GAP_PT * (len(pieces) - 1)
    x_pt = FLOW_CENTER_X_PT - total / 2
    for piece, width_pt in zip(pieces, widths, strict=True):
        artist = container.text(
            x_pt,
            y_pt,
            piece.text,
            ha="left",
            va="center",
            fontsize=EQUATION_FONT_PT,
            color=piece.color,
            zorder=5,
        )
        if piece.chip is not None:
            draw_chip(container, artist, piece.chip)
        x_pt += width_pt + EQUATION_LINE_GAP_PT


def draw_chip(container: Axes, artist: Text, chip: tuple[float, float, float, float]) -> None:
    """Draw a rounded tint behind one parameter name in the equation."""
    renderer = cast(FigureCanvasAgg, container.figure.canvas).get_renderer()
    bbox = artist.get_window_extent(renderer)
    inverse = container.transData.inverted()
    (x0, y0), (x1, y1) = inverse.transform([[bbox.x0, bbox.y0], [bbox.x1, bbox.y1]])
    container.add_patch(
        FancyBboxPatch(
            (x0 - CHIP_PAD_PT, y0 - CHIP_PAD_PT),
            (x1 - x0) + 2 * CHIP_PAD_PT,
            (y1 - y0) + 2 * CHIP_PAD_PT,
            boxstyle=f"round,pad=0,rounding_size={CHIP_ROUNDING_PT}",
            facecolor=chip,
            edgecolor="none",
            zorder=4,
        )
    )


def draw_parameter_names(container: Axes) -> None:
    """Draw the extracted parameters in their colours, centred in their box."""
    widths = [
        text_extent_pt(container, text, PARAMETER_NAME_FONT_PT)[0] for text, _ in PARAMETER_NAMES
    ]
    total = sum(widths) + PARAMETER_NAME_GAP_PT * (len(PARAMETER_NAMES) - 1)
    x_pt = FLOW_CENTER_X_PT - total / 2
    for (text, color), width_pt in zip(PARAMETER_NAMES, widths, strict=True):
        container.text(
            x_pt,
            PARAMETER_BOX_Y_PT + 6.5,
            text,
            ha="left",
            va="center",
            fontsize=PARAMETER_NAME_FONT_PT,
            color=color,
            fontweight="bold",
        )
        x_pt += width_pt + PARAMETER_NAME_GAP_PT


def fit_curve(fit: GeneFit) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the fitted curve on a dense grid."""
    curve_x = np.linspace(fit.generations.min(), fit.generations.max(), CURVE_SAMPLES)
    return curve_x, gompertz(curve_x, fit.amplitude, fit.rate, fit.lag)


def panel_x_limits(fit: GeneFit) -> tuple[float, float]:
    """Return the shared x limits around the generation grid."""
    span = np.ptp(fit.generations)
    return fit.generations.min() - X_PADDING * span, fit.generations.max() + X_PADDING * span


def panel_b_x_limits(fit: GeneFit) -> tuple[float, float]:
    """Return panel B's x limits, with room left of generation 0 for a DL = 0 guide."""
    span = np.ptp(fit.generations)
    return (
        fit.generations.min() - B_X_PADDING_LEFT * span,
        fit.generations.max() + B_X_PADDING_RIGHT * span,
    )


def draw_fit_series(ax: Axes, fit: GeneFit, curve_x: np.ndarray, curve_y: np.ndarray) -> None:
    """Draw the fitted curve, the measured points and the zero baseline."""
    ax.plot(curve_x, curve_y, color=FIT_COLOR, linewidth=1.5, zorder=3)
    ax.scatter(
        fit.generations,
        fit.observed,
        s=OBSERVED_MARKER_SIZE,
        color=POINT_COLOR,
        edgecolors="white",
        linewidths=OBSERVED_MARKER_EDGE_PT,
        zorder=4,
    )
    ax.axhline(0, color=FURNITURE_COLOR, linestyle=(0, (1, 2)), linewidth=ZERO_LINE_WIDTH_PT)


def render_anatomy_plot(ax: Axes, fit: GeneFit) -> None:
    """Draw the example fit with A, DR and DL written onto their geometry."""
    curve_x, curve_y = fit_curve(fit)
    low, high = y_data_bounds(fit, curve_y)
    y_bottom, y_top = padded_y_limits(low, high)
    x_left, x_right = panel_x_limits(fit)

    draw_fit_series(ax, fit, curve_x, curve_y)

    ax.set_xlim(x_left, x_right)
    ax.set_ylim(y_bottom, y_top)
    draw_plateau_guide(ax, fit, x_left, x_right)
    draw_rate_guide(ax, fit, x_left, x_right, y_top, y_bottom)
    draw_delay_guide(ax, fit, y_bottom)
    ax.set_xticks(fit.generations, display_ticks(fit.generations))
    ax.set_yticks(y_tick_values(low, high))
    ax.tick_params(labelsize=PLOT_TICK_FONT_PT)
    ax.set_xlabel(X_LABEL, fontsize=PLOT_LABEL_FONT_PT)
    ax.set_ylabel(Y_LABEL, fontsize=PLOT_LABEL_FONT_PT)


def draw_plateau_guide(
    ax: Axes,
    fit: GeneFit,
    x_left: float,
    x_right: float,
    *,
    with_label: bool = True,
    dot_size: float = GUIDE_DOT_SIZE,
    line_width: float = GUIDE_LINE_WIDTH_PT,
) -> None:
    """Draw the plateau line, its end dot and its label."""
    ax.axhline(fit.amplitude, color=AMPLITUDE_COLOR, linestyle=(0, (4, 2)), linewidth=line_width)
    ax.scatter([x_right - 0.3], [fit.amplitude], s=dot_size, color=AMPLITUDE_COLOR, zorder=5)
    if not with_label:
        return
    # Above the line on the left: the tangent and the late-time curve occupy the
    # line's right half, and the delay line sits far enough right of this label.
    ax.text(
        x_left + 0.3,
        fit.amplitude + 0.12,
        PLATEAU_LABEL,
        color=AMPLITUDE_COLOR,
        ha="left",
        va="bottom",
        fontsize=PLOT_TICK_FONT_PT,
    )


def draw_rate_guide(
    ax: Axes,
    fit: GeneFit,
    x_left: float,
    x_right: float,
    y_top: float,
    y_bottom: float,
    *,
    with_label: bool = True,
    dot_size: float = GUIDE_DOT_SIZE,
    line_width: float = GUIDE_LINE_WIDTH_PT,
) -> None:
    """Draw the inflection tangent, its inflection dot and the slope label."""
    tangent_x = _tangent_segment(fit, x_left, x_right, y_top, y_bottom)
    ax.plot(
        tangent_x,
        fit.rate * (tangent_x - fit.lag),
        color=RATE_COLOR,
        linestyle=(0, (5, 1.4, 0.8, 1.4)),
        linewidth=line_width,
        zorder=2,
    )
    ax.scatter([fit.t_inflection], [fit.y_inflection], s=dot_size, color=RATE_COLOR, zorder=5)
    if not with_label:
        return
    ax.text(
        3.3,
        1.15,
        RATE_LABEL,
        color=RATE_COLOR,
        ha="left",
        va="bottom",
        fontsize=PLOT_TICK_FONT_PT,
    )


def _tangent_segment(
    fit: GeneFit, x_left: float, x_right: float, y_top: float, y_bottom: float
) -> np.ndarray:
    """Clip the tangent line to the axes box, returning its two x endpoints."""
    crossing_top = fit.lag + y_top / fit.rate
    crossing_bottom = fit.lag + y_bottom / fit.rate
    start = max(min(crossing_top, crossing_bottom), x_left)
    end = min(max(crossing_top, crossing_bottom), x_right)
    return np.array([start, end])


def draw_delay_guide(
    ax: Axes,
    fit: GeneFit,
    y_bottom: float,
    *,
    with_label: bool = True,
    dot_size: float = GUIDE_DOT_SIZE,
    line_width: float = GUIDE_LINE_WIDTH_PT,
) -> None:
    """Draw the delay line, its two dots and the intercept label."""
    line_bottom = y_bottom + DELAY_BOTTOM_MARGIN
    ax.plot(
        [fit.lag, fit.lag],
        [0, line_bottom],
        color=LAG_COLOR,
        linestyle=(0, (1, 1.6)),
        linewidth=line_width,
        zorder=2,
    )
    ax.scatter([fit.lag, fit.lag], [0, line_bottom], s=dot_size, color=LAG_COLOR, zorder=5)
    if not with_label:
        return
    ax.text(
        fit.lag + 0.25,
        0.12,
        LAG_LABEL,
        color=LAG_COLOR,
        ha="left",
        va="bottom",
        fontsize=PLOT_TICK_FONT_PT,
    )


# =============================================================================
# 6. ASSEMBLY
# =============================================================================
def mount_on_a4(canvas_pdf: Path, output_pdf: Path) -> None:
    """Inset the content canvas from the A4 page's top-left corner; the rest stays blank."""
    canvas_page = PdfReader(str(canvas_pdf)).pages[0]
    content_height_pt = float(canvas_page.mediabox.height)
    bottom_pt = A4_HEIGHT_PT - PAGE_MARGIN_PT - content_height_pt
    writer = PdfWriter()
    page = writer.add_blank_page(width=A4_WIDTH_PT, height=A4_HEIGHT_PT)
    page.merge_transformed_page(canvas_page, Transformation().translate(PAGE_MARGIN_PT, bottom_pt))
    writer.write(str(output_pdf))


def rasterize_review(pdf: Path, review_png: Path) -> None:
    """Rasterize the assembled PDF at the house review resolution."""
    subprocess.run(
        [
            "pdftocairo",
            "-png",
            "-r",
            str(HOUSE_FIGURE_DPI),
            "-singlefile",
            str(pdf),
            str(review_png.with_suffix("")),
        ],
        check=True,
    )


# =============================================================================
# 7. IMPERATIVE SHELL (boundary validation, orchestration, CLI)
# =============================================================================
def read_tsv_columns(path: Path, required: tuple[str, ...]) -> pd.DataFrame:
    """Read only the required TSV columns and fail once if any is missing."""
    table = pd.read_csv(path, sep="\t", usecols=lambda column: column in required)
    missing_columns = set(required) - set(table.columns)
    if missing_columns:
        raise ValueError(f"{path} is missing required columns: {sorted(missing_columns)}")
    return table


def run(config: Figure2Config) -> tuple[float, float]:
    """Build Figure 2 end to end and return the page's (width, height) in points."""
    config.validate()
    apply_house_style()

    fitting = read_tsv_columns(config.gene_fitting, GENE_FIT_COLUMNS)
    distributions = build_parameter_distributions(
        read_tsv_columns(config.gene_fitting, DISTRIBUTION_COLUMNS)
    )
    fits = {gene: build_gene_fit(fitting, gene) for gene in (PANEL_A_GENE, *PANEL_B_GENES)}

    mp = draw_figure(fits[PANEL_A_GENE], tuple(fits[gene] for gene in PANEL_B_GENES), distributions)
    assert mp.fig is not None
    # cnsplots sets savefig.bbox to "tight"; pass bbox_inches=None so the page
    # keeps the multipanel canvas instead of cropping to the drawn artists.
    cns.savefig(config.canvas_pdf, fig=mp.fig, bbox_inches=None)
    mount_on_a4(config.canvas_pdf, config.output_pdf)
    rasterize_review(config.output_pdf, config.output_review)

    logger.success(
        f"Wrote {config.output_pdf} (A4 {A4_WIDTH_PT:.2f} x {A4_HEIGHT_PT:.2f} pt) "
        f"and {config.output_review}"
    )
    return A4_WIDTH_PT, A4_HEIGHT_PT


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Build the composite Figure 2")
    parser.add_argument(
        "--gene-fitting", type=Path, required=True, help="Release gene-level fitting TSV"
    )
    parser.add_argument("--output-pdf", type=Path, required=True, help="Final figure PDF")
    parser.add_argument("--output-review", type=Path, required=True, help="Review PNG")
    parser.add_argument("--canvas", type=Path, required=True, help="Canvas PDF before the A4 mount")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging"
    )
    return parser.parse_args()


def setup_logger(verbose: bool) -> None:
    """Configure the loguru sink at INFO or DEBUG."""
    logger.remove()
    logger.add(sys.stderr, level="DEBUG" if verbose else "INFO")


def main() -> int:
    """Build the config, render the figure, and report the outcome."""
    args = parse_args()
    setup_logger(args.verbose)
    config = Figure2Config(
        gene_fitting=args.gene_fitting,
        output_pdf=args.output_pdf,
        output_review=args.output_review,
        canvas_pdf=args.canvas,
    )
    try:
        run(config)
    except Exception as e:
        logger.exception(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
