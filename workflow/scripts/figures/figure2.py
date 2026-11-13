#!/usr/bin/env python3

"""
Figure 2 — Depletion kinetics: the Gompertz fit behind DR and DL
=================================================================

Panel A is the worked example that introduces the curve model every downstream
DR/DL feature is read from: one essential gene (nrs1), its five measured
gene-level LFCs through the 12.48-generation time course, and the release's
Gompertz fit. The fitted curve is drawn from the released amplitude, rate and
lag parameters (A, DR, DL) at a resolution fine enough to look smooth; the
release's five ``*_fitted`` columns are used only to confirm that this script's
model reproduces the fit of record, never to place the curve. The formula and
the four fitted values are written inside the panel.

The gene is fixed by name (``GENE_NAME``): it is the example the figure tells
its story with, not a per-run choice. Every drawing decision lives in this
file; the script imports no stage module from ``workflow/src`` and calls no
upstream renderer, so retuning Figure 2 cannot change another stage's output.

The page is an A4 manuscript page (210 x 297 mm). ``cns.multipanel`` sizes only
the canvas width and lets the rows decide the height, so the flow below draws a
content-sized canvas and the assembly mounts it 40 pt in from the A4 page's
top-left corner, leaving the rest of the page blank. Each panel's box is its
axes size plus the y-axis decoration multipanel measures at render time, so a
wide ylabel pushes only its own panel right; the letters are placed by
multipanel at the left edge of each panel's column.

Input
-----
- <HD_DIT_HAP release>/gene_level/gene_level_fitting_statistics.tsv
  Columns: Name, time_points, A, DR, DL, R2, YES0..YES4, YES0_fitted..YES4_fitted.

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
Version:  1.0.0
"""

# =============================================================================
# 1. IMPORTS (Standard Library -> Third-Party)
# =============================================================================
import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
import sys

import cnsplots as cns
from loguru import logger
from matplotlib.axes import Axes
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
BELOW_DECORATION_PT = 26.0  # x ticks + xlabel
# 120 pt is the example panel's width: wide enough for the one-line formula at
# the house 7 pt without shrinking the type below it.
PANEL_WIDTH_PT = 120.0
PANEL_HEIGHT_PT = 100.0

# --- Panel A vocabulary -----------------------------------------------------
GENE_NAME = "nrs1"
NAME_COLUMN = "Name"
TIME_POINTS_COLUMN = "time_points"
AMPLITUDE_COLUMN = "A"
RATE_COLUMN = "DR"
LAG_COLUMN = "DL"
R2_COLUMN = "R2"
VALUE_COLUMNS = ("YES0", "YES1", "YES2", "YES3", "YES4")
FITTED_COLUMNS = tuple(f"{column}_fitted" for column in VALUE_COLUMNS)
# The release keeps only two or three decimals for A/DR/DL and the fitted
# values, so evaluating the model from the stored parameters can sit a few
# units in the last place away from the stored fit; 5e-3 covers the rounding
# of every parameter without accepting a different model.
FITTED_TOLERANCE = 5e-3

# Observed points take the accent, the fitted curve the greyscale-separable
# teal of slot 1: shape and colour both separate the two series.
OBSERVED_COLOR = cns.get_hexcolors_from_apalette([0], HOUSE_PALETTE)[0]
FIT_COLOR = cns.get_hexcolors_from_apalette([1], HOUSE_PALETTE)[0]
OBSERVED_MARKER_SIZE = 24.0
OBSERVED_MARKER_EDGE_PT = 0.6
FIT_LINE_WIDTH_PT = 1.5
CURVE_SAMPLES = 300

X_LABEL = "Generations"
Y_LABEL = r"Gene-level LFC ($\log_2$)"
# One decimal is the display precision the release's generation grid needs
# (2.352 -> 2.4); a trailing ".0" is dropped so the first tick reads 0.
TICK_DECIMALS = 1

# The data fill the lower ~70% of the axes; the top band is left for the model
# block. Padding is a fraction of the combined observed/fitted span so the
# panel keeps its proportions if the example gene ever changes.
Y_PADDING_BOTTOM = 0.08
Y_PADDING_TOP = 0.30
X_PADDING = 0.03
# Two LFC units per tick across the data range only: the reserved annotation
# band above the data carries no tick, so it reads as headroom, not as scale.
Y_TICK_STEP = 2.0

# --- Panel A annotation -----------------------------------------------------
# The exact model as implemented upstream: alpha = DR * e / A. Written as one
# display equation; the parameters below carry the fitted values.
FIT_FORMULA = r"$y=A\,\exp\left(-\exp\left[\dfrac{DR\cdot e}{A}\left(DL-x\right)+1\right]\right)$"
ANNOTATION_X = 0.985  # axes fraction, right-aligned
ANNOTATION_TOP = 0.975
ANNOTATION_LINE_STEP = 0.145  # axes fraction between the three lines
ANNOTATION_FONT_PT = 7.0

# --- Input contract ---------------------------------------------------------
GENE_FIT_COLUMNS = (
    NAME_COLUMN,
    TIME_POINTS_COLUMN,
    AMPLITUDE_COLUMN,
    RATE_COLUMN,
    LAG_COLUMN,
    R2_COLUMN,
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


def build_gene_fit(fitting: pd.DataFrame) -> GeneFit:
    """Select the example gene and tie this script's model to the released fit."""
    is_example = fitting[NAME_COLUMN] == GENE_NAME
    matches = fitting.loc[is_example]
    if matches.empty:
        raise ValueError(f"{GENE_NAME} is not in the gene-level fitting table")
    if len(matches) > 1:
        raise ValueError(f"{GENE_NAME} appears {len(matches)} times in the fitting table")

    row = matches.iloc[0]
    generations = np.array(
        [float(point) for point in row[TIME_POINTS_COLUMN].split(",")], dtype=float
    )
    if len(generations) != len(VALUE_COLUMNS):
        raise ValueError(
            f"{GENE_NAME} has {len(generations)} time points for {len(VALUE_COLUMNS)} values"
        )

    amplitude = float(row[AMPLITUDE_COLUMN])
    rate = float(row[RATE_COLUMN])
    lag = float(row[LAG_COLUMN])
    released_fitted = row[list(FITTED_COLUMNS)].to_numpy(dtype=float)
    model_fitted = gompertz(generations, amplitude, rate, lag)
    deviation = float(np.max(np.abs(model_fitted - released_fitted)))
    if deviation > FITTED_TOLERANCE:
        raise ValueError(
            f"{GENE_NAME}: model from A/DR/DL deviates {deviation:.4f} from the released fit"
        )

    return GeneFit(
        name=str(row[NAME_COLUMN]),
        generations=generations,
        observed=row[list(VALUE_COLUMNS)].to_numpy(dtype=float),
        amplitude=amplitude,
        rate=rate,
        lag=lag,
        r2=float(row[R2_COLUMN]),
    )


def display_ticks(generations: np.ndarray) -> list[str]:
    """Format the generation grid at the panel's display precision."""
    labels = []
    for generation in generations:
        label = f"{generation:.{TICK_DECIMALS}f}"
        labels.append(label.removesuffix(".0"))
    return labels


def y_data_bounds(fit: GeneFit, curve: np.ndarray) -> tuple[float, float]:
    """Return the lowest and highest y among the points and the drawn curve."""
    low = min(float(fit.observed.min()), float(curve.min()))
    high = max(float(fit.observed.max()), float(curve.max()))
    return low, high


def padded_y_limits(low: float, high: float) -> tuple[float, float]:
    """Pad the data bounds, leaving the wider band on top for the model block."""
    span = high - low
    return low - Y_PADDING_BOTTOM * span, high + Y_PADDING_TOP * span


def y_tick_values(low: float, high: float) -> np.ndarray:
    """Tick the 2-unit grid that falls inside the data bounds."""
    first = np.ceil(low / Y_TICK_STEP) * Y_TICK_STEP
    last = np.floor(high / Y_TICK_STEP) * Y_TICK_STEP
    return np.arange(first, last + Y_TICK_STEP / 2, Y_TICK_STEP)


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


def draw_figure(fit: GeneFit) -> cns.multipanel:
    """Draw panel A and return the multipanel that owns the page.

    Returns the multipanel itself, not just its figure: cnsplots keeps the
    measured-layout draw callback as a weak reference, so the caller must hold
    the object until the file is written or the panels stay unmeasured.
    """
    mp = cns.multipanel(max_width=CANVAS_WIDTH_PT, unit="pt")
    ax = mp.panel(
        "A",
        width=PANEL_WIDTH_PT,
        height=PANEL_HEIGHT_PT,
        margin_bottom=BELOW_DECORATION_PT,
    )
    render_gene_fit_panel(ax, fit)
    return mp


def render_gene_fit_panel(ax: Axes, fit: GeneFit) -> None:
    """Draw the measured points, the smooth fitted curve and the model block."""
    curve_x = np.linspace(fit.generations.min(), fit.generations.max(), CURVE_SAMPLES)
    curve_y = gompertz(curve_x, fit.amplitude, fit.rate, fit.lag)

    ax.scatter(
        fit.generations,
        fit.observed,
        s=OBSERVED_MARKER_SIZE,
        color=OBSERVED_COLOR,
        edgecolors="white",
        linewidths=OBSERVED_MARKER_EDGE_PT,
        zorder=3,
    )
    ax.plot(curve_x, curve_y, color=FIT_COLOR, linewidth=FIT_LINE_WIDTH_PT, zorder=2)

    low, high = y_data_bounds(fit, curve_y)
    ax.set_xlim(
        fit.generations.min() - X_PADDING * np.ptp(fit.generations),
        fit.generations.max() + X_PADDING * np.ptp(fit.generations),
    )
    ax.set_ylim(*padded_y_limits(low, high))
    ax.set_xticks(fit.generations, display_ticks(fit.generations))
    ax.set_yticks(y_tick_values(low, high))
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)
    ax.set_title(fit.name, fontstyle="italic", fontweight="bold")
    add_model_block(ax, fit)


def add_model_block(ax: Axes, fit: GeneFit) -> None:
    """Write the fitted formula and values in the axes' free top-right band."""
    # A zero lag is stored as -0.0 by the release's three-decimal round trip;
    # print the value, not its signed zero.
    lag = 0.0 if fit.lag == 0 else fit.lag
    lines = (
        FIT_FORMULA,
        rf"$A={fit.amplitude:.2f},\quad DR={fit.rate:.2f}$",
        rf"$DL={lag:.2f},\quad R^2={fit.r2:.3f}$",
    )
    for index, line in enumerate(lines):
        ax.text(
            ANNOTATION_X,
            ANNOTATION_TOP - index * ANNOTATION_LINE_STEP,
            line,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=ANNOTATION_FONT_PT,
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
    fit = build_gene_fit(fitting)

    mp = draw_figure(fit)
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
