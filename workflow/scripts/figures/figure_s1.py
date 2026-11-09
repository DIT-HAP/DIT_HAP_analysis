#!/usr/bin/env python3

"""
Figure S1 — Library-prep reproducibility: PBL/PBR, technical and biological replicates
======================================================================================

Three density scatter panels from the 1a_pcr_qc stage, plus the distribution of
insertions per gene for HD_DIT_HAP:

  A  PBL versus PBR reads of one library (the two primers agree);
  B  technical replicate — one sample processed in two upstream projects;
  C  biological replicate — two samples within one project;
  D  how many insertions each coding gene carries (histogram of counts).

The 1a_pcr_qc stage first drew panels A-C inside its 2x2 QC figure. This script
re-renders them on a publication canvas: it imports no stage module from
``workflow/src`` and carries its own copy of the density colouring, log axes
and n/r/P annotation, so retuning Figure S1 cannot change another stage's
output. Panel D is new here: it counts in-gene insertions per gene straight
from the released insertion annotations. Genes with no in-gene insertion are
absent from that table, so the distribution covers only the genes that carry at
least one.

Input
-----
- results/1a_pcr_qc/_work/pbl_pbr.parquet   Columns: PBL, PBR
- results/1a_pcr_qc/_work/tech.parquet      Columns: Reads_1, Reads_2
- results/1a_pcr_qc/_work/bio.parquet       Columns: Reads_1, Reads_2
- <HD_DIT_HAP release>/insertion_level/annotations.tsv.gz
  Columns: Type, Systematic ID (panel D)

Output
------
- results/7a_figures/figure_s1.pdf
- results/7a_figures/figure_s1.review.png

Usage
-----
    python workflow/scripts/figures/figure_s1.py \\
        --pbl-pbr results/1a_pcr_qc/_work/pbl_pbr.parquet \\
        --tech results/1a_pcr_qc/_work/tech.parquet \\
        --bio results/1a_pcr_qc/_work/bio.parquet \\
        --annotations <release>/insertion_level/annotations.tsv.gz \\
        --output-pdf results/7a_figures/figure_s1.pdf \\
        --output-review results/7a_figures/figure_s1.review.png

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-09
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
from matplotlib.collections import PathCollection
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import gaussian_kde, pearsonr

# =============================================================================
# 2. CONSTANTS (house style, layout, draw vocabulary)
# =============================================================================
# --- House style (self-contained copy; see docs/figures.md) -----------------
HOUSE_PALETTE = "Cell"
DENSITY_CMAP = "viridis"
FONT_SANS_SERIF = (
    "Arial",
    "Helvetica",
    "Nimbus Sans",
    "Liberation Sans",
    "DejaVu Sans",
)
HOUSE_FIGURE_DPI = 360
FURNITURE_COLOR = cns.GRAY  # guides and annotations, never data
# Panel D has one bar series; palette position 2 is the greyscale-safe series slot.
COUNT_BAR_COLOR = cns.get_hexcolors_from_apalette([2], HOUSE_PALETTE)[0]

# --- Canvas and panel geometry, in points -----------------------------------
CANVAS_WIDTH_PT = 493.0
BOTTOM_MARGIN_PT = 2.0
PANEL_LETTER_ROOM_PT = 10.0  # clear band above the row for its letters
PANEL_LETTER_GAP_PT = 2.0
PANEL_LETTER_FONT_PT = 8.0
BELOW_DECORATION_PT = 30.0  # x ticks + xlabel

# Three square panels across the page. 120 pt sits above the house 100 pt
# legible floor, and the pitch leaves each panel its own y tick labels and
# ylabel.
PANEL_SIZE_PT = 120.0
PANEL_PITCH_PT = 160.0
PANEL_FIRST_X_PT = 30.0

# Panel D sits on a second row at the house SQUARE footprint, the same 100 pt
# side as the panel contract's baseline.
ROW_GAP_PT = 10.0
PANEL_D_X_PT = 30.0
PANEL_D_WIDTH_PT = 100.0
PANEL_D_HEIGHT_PT = 100.0

# --- Density scatter --------------------------------------------------------
# Point size and alpha are the dense-cloud exception from the house style:
# cnsplots' s=7 fills in once thousands of points overlap.
DENSITY_SCATTER_S = 1.5
DENSITY_SCATTER_ALPHA = 0.6
DENSITY_CBAR_BOUNDS = (0.82, 0.02, 0.035, 0.30)
DENSITY_CBAR_LABEL = "Density"

# KDE cost: gaussian_kde is O(n_eval x n_train), so above DENSITY_EXACT_MAX_N
# the KDE is fitted on a fixed-seed subsample, evaluated on a grid and
# interpolated back. 5000 points is ~2.5e7 kernel evaluations.
DENSITY_EXACT_MAX_N = 5000
DENSITY_KDE_SUBSAMPLE_N = 5000
DENSITY_GRID_SIZE = 128
DENSITY_SEED = 0  # fixed so repeat renders are identical

# --- Statistics annotation --------------------------------------------------
STATS_X = 0.05
STATS_Y = 0.95
P_UNDERFLOW_EXPONENT = -9  # below 1E-9 the text reports "P < 1E-9" instead of an exponent
P_UNDERFLOW = 10.0**P_UNDERFLOW_EXPONENT
P_UNDERFLOW_TEXT = rf"$P<1\mathrm{{E}}{{{P_UNDERFLOW_EXPONENT}}}$"
MINOR_TICK_LENGTH_RATIO = 0.57  # minor vs major tick length

# --- Panel D vocabulary -----------------------------------------------------
# Counts are taken where the release annotation says the insertion lands in a
# coding gene; a gene absent from that table has no in-gene insertion and is
# therefore outside the distribution. The linear bins stop at COUNT_TAIL_START
# and everything above shares one bar, labelled so the last bar is not read as
# another bin-width.
ANNOTATION_TYPE_COLUMN = "Type"
ANNOTATION_GENE_COLUMN = "Systematic ID"
ANNOTATION_COLUMNS = (ANNOTATION_TYPE_COLUMN, ANNOTATION_GENE_COLUMN)
GENE_TYPE = "Coding gene"
COUNT_BIN_WIDTH = 25
COUNT_TAIL_START = 300  # the last linear bin edge
COUNT_TAIL_LABEL = f"≥{COUNT_TAIL_START}"
COUNT_LABEL_EVERY = 3  # 13 bin edges do not fit under a 100 pt axes
COUNT_X_LABEL = "Insertions per gene"
COUNT_Y_LABEL = "Number of genes"
COUNT_PANEL_LETTER = "D"


# --- Panel vocabulary -------------------------------------------------------
@dataclass(kw_only=True, frozen=True)
class ScatterSpec:
    """One panel: which columns to compare and how to label the axes."""

    x: str
    y: str
    xlabel: str
    ylabel: str


PANEL_SPECS = (
    ScatterSpec(x="PBL", y="PBR", xlabel="PBL Reads", ylabel="PBR Reads"),
    ScatterSpec(
        x="Reads_1",
        y="Reads_2",
        xlabel="Reads of Technical Replicate 1",
        ylabel="Reads of Technical Replicate 2",
    ),
    ScatterSpec(
        x="Reads_1",
        y="Reads_2",
        xlabel="Reads of Biological Replicate 1",
        ylabel="Reads of Biological Replicate 2",
    ),
)
PANEL_LETTERS = ("A", "B", "C")


# =============================================================================
# 3. CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class FigureS1Config:
    """Resolved input/output paths for Figure S1."""

    pbl_pbr: Path
    tech: Path
    bio: Path
    annotations: Path
    output_pdf: Path
    output_review: Path

    def validate(self) -> None:
        """Raise if a required input is missing, then create the output dirs."""
        inputs = [self.pbl_pbr, self.tech, self.bio, self.annotations]
        for path in inputs:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for output in [self.output_pdf, self.output_review]:
            output.parent.mkdir(parents=True, exist_ok=True)


@dataclass(kw_only=True, frozen=True)
class Layout:
    """Resolved point geometry of the two panel rows."""

    canvas_height_pt: float
    panel_xs_pt: tuple[float, ...]
    panel_rects: tuple[tuple[float, float, float, float], ...]
    panel_d_top_pt: float
    panel_d_rect: tuple[float, float, float, float]


# =============================================================================
# 4. DATA PREPARATION (pure; schemas are validated at the I/O boundary)
# =============================================================================
def build_scatter_frame(table: pd.DataFrame, spec: ScatterSpec) -> pd.DataFrame:
    """Keep the pairs both axes can draw: strictly positive read counts."""
    is_positive = (table[spec.x] > 0) & (table[spec.y] > 0)
    return table.loc[is_positive, [spec.x, spec.y]]


def build_insertion_counts(annotations: pd.DataFrame) -> pd.Series:
    """Count in-gene insertions per gene."""
    is_in_gene = annotations[ANNOTATION_TYPE_COLUMN] == GENE_TYPE
    return annotations.loc[is_in_gene, ANNOTATION_GENE_COLUMN].value_counts()


def build_count_histogram(counts: pd.Series) -> pd.DataFrame:
    """Bin per-gene insertion counts; the long tail shares one labelled bar."""
    linear_edges = range(0, COUNT_TAIL_START, COUNT_BIN_WIDTH)
    labels = [str(edge) for edge in linear_edges] + [COUNT_TAIL_LABEL]
    edges = [*linear_edges, COUNT_TAIL_START, np.inf]
    binned = pd.cut(counts, bins=edges, right=False, labels=labels)
    return binned.value_counts(sort=False).rename_axis("bin_label").reset_index(name="n_genes")


def build_layout() -> Layout:
    """Derive the canvas height and every panel box from the two-row geometry."""
    row_a_bottom = PANEL_LETTER_ROOM_PT + PANEL_SIZE_PT + BELOW_DECORATION_PT
    panel_d_top = row_a_bottom + ROW_GAP_PT + PANEL_LETTER_ROOM_PT
    canvas_height = panel_d_top + PANEL_D_HEIGHT_PT + BELOW_DECORATION_PT + BOTTOM_MARGIN_PT

    panel_xs = tuple(PANEL_FIRST_X_PT + index * PANEL_PITCH_PT for index in range(len(PANEL_SPECS)))
    panel_rects = tuple(
        _axes_rect(x_pt, PANEL_LETTER_ROOM_PT, PANEL_SIZE_PT, PANEL_SIZE_PT, canvas_height)
        for x_pt in panel_xs
    )
    return Layout(
        canvas_height_pt=canvas_height,
        panel_xs_pt=panel_xs,
        panel_rects=panel_rects,
        panel_d_top_pt=panel_d_top,
        panel_d_rect=_axes_rect(
            PANEL_D_X_PT, panel_d_top, PANEL_D_WIDTH_PT, PANEL_D_HEIGHT_PT, canvas_height
        ),
    )


def _axes_rect(
    x_pt: float, y_top_pt: float, width_pt: float, height_pt: float, canvas_height_pt: float
) -> tuple[float, float, float, float]:
    """Convert a top-left point box into matplotlib figure fractions."""
    bottom_pt = canvas_height_pt - y_top_pt - height_pt
    return (
        x_pt / CANVAS_WIDTH_PT,
        bottom_pt / canvas_height_pt,
        width_pt / CANVAS_WIDTH_PT,
        height_pt / canvas_height_pt,
    )


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
    cns.settings.palette_seq = DENSITY_CMAP
    cns.setup_matplotlib(color_cycle=HOUSE_PALETTE, color_map=DENSITY_CMAP)


def draw_figure(
    frames: tuple[pd.DataFrame, ...], count_histogram: pd.DataFrame, layout: Layout
) -> Figure:
    """Draw the three scatter panels, panel D, and every panel letter."""
    fig = cns.figure(width=CANVAS_WIDTH_PT, height=layout.canvas_height_pt, unit="pt")
    panels = zip(PANEL_LETTERS, PANEL_SPECS, frames, strict=True)

    for index, (letter, spec, frame) in enumerate(panels):
        ax = fig.add_axes(layout.panel_rects[index])
        render_replicate_panel(ax, frame, spec)
        add_panel_letter(
            fig, letter, layout.panel_xs_pt[index], PANEL_LETTER_ROOM_PT, layout.canvas_height_pt
        )

    ax_d = fig.add_axes(layout.panel_d_rect)
    render_insertion_count_panel(ax_d, count_histogram)
    add_panel_letter(
        fig, COUNT_PANEL_LETTER, PANEL_D_X_PT, layout.panel_d_top_pt, layout.canvas_height_pt
    )

    return fig


def add_panel_letter(
    fig: Figure, label: str, x_pt: float, y_top_pt: float, canvas_height_pt: float
) -> None:
    """Place a bold panel letter just above the panel's top-left corner."""
    fig.text(
        x_pt / CANVAS_WIDTH_PT,
        1.0 - (y_top_pt - PANEL_LETTER_GAP_PT) / canvas_height_pt,
        label,
        ha="left",
        va="bottom",
        fontweight="bold",
        fontsize=PANEL_LETTER_FONT_PT,
    )


def render_replicate_panel(ax: Axes, frame: pd.DataFrame, spec: ScatterSpec) -> None:
    """Draw one log-log density scatter with its identity guide and n/r/P block."""
    x = frame[spec.x].to_numpy()
    y = frame[spec.y].to_numpy()
    log_x = np.log10(x)
    log_y = np.log10(y)

    draw_density_scatter(ax, x, y, point_density(log_x, log_y))
    draw_identity_line(ax, x, y)
    annotate_fit_stats(ax, log_x, log_y)

    ax.set_xlabel(spec.xlabel)
    ax.set_ylabel(spec.ylabel)
    apply_log_scale(ax)


def render_insertion_count_panel(ax: Axes, frame: pd.DataFrame) -> None:
    """Draw the per-gene insertion count histogram, one bar per count bin."""
    cns.barplot(data=frame, x="bin_label", y="n_genes", color=COUNT_BAR_COLOR, ax=ax)
    # Keep every other bin edge, last one included, so the labels cannot collide.
    ax.set_xticks(ax.get_xticks()[::COUNT_LABEL_EVERY])
    ax.set_xlabel(COUNT_X_LABEL)
    ax.set_ylabel(COUNT_Y_LABEL)


def draw_density_scatter(ax: Axes, x: np.ndarray, y: np.ndarray, density: np.ndarray) -> None:
    """Colour the points by 2D density and inset the low/high colorbar."""
    # Ascending, so the crowded points draw last and are not buried under the
    # sparse ones.
    order = np.argsort(density)
    ranked = density[order]
    span = np.ptp(ranked)
    colour = (ranked - ranked.min()) / span

    # [0, 1]-normalized so the colorbar's ticks sit at the bar ends; density is
    # renormalised per panel, so only the ordering carries meaning.
    mappable = ax.scatter(
        x[order],
        y[order],
        c=colour,
        norm=Normalize(vmin=0, vmax=1),
        s=DENSITY_SCATTER_S,
        alpha=DENSITY_SCATTER_ALPHA,
        linewidths=0,
        rasterized=True,
    )
    add_density_colorbar(ax, mappable)


def add_density_colorbar(ax: Axes, mappable: PathCollection) -> None:
    """Inset a thin vertical colorbar labelled low/high inside the panel's bottom-right corner."""
    cax = ax.inset_axes(DENSITY_CBAR_BOUNDS)
    colorbar = ax.figure.colorbar(mappable, cax=cax, ticks=[0, 1])

    labels = colorbar.ax.set_yticklabels(["low", "high"])
    # End ticks anchor their label on the tick centre, pushing half the text
    # past the bar; aligning each label inwards keeps both inside the panel.
    labels[0].set_verticalalignment("bottom")
    labels[-1].set_verticalalignment("top")
    # Above the bar, not beside it: the bar is far narrower than the word.
    colorbar.ax.set_title(DENSITY_CBAR_LABEL)


def draw_identity_line(ax: Axes, x: np.ndarray, y: np.ndarray) -> None:
    """Draw the dashed x = y guide across the data range."""
    low = min(x.min(), y.min())
    high = max(x.max(), y.max())
    ax.plot([low, high], [low, high], color=FURNITURE_COLOR, linestyle="--")


def annotate_fit_stats(ax: Axes, log_x: np.ndarray, log_y: np.ndarray) -> None:
    """Write n, Pearson r and its P value in the panel's top-left corner."""
    r, p = pearsonr(log_x, log_y)
    p_text = P_UNDERFLOW_TEXT if p < P_UNDERFLOW else _p_value_text(p)
    text = f"$n={log_x.size}$\n$r={r:.2f}$\n{p_text}"
    ax.text(STATS_X, STATS_Y, text, transform=ax.transAxes, ha="left", va="top")


def _p_value_text(p: float) -> str:
    """Format P as mantissa E exponent, e.g. 3.80E-29."""
    exponent = int(np.floor(np.log10(p)))
    mantissa = p / 10**exponent
    return rf"$P={mantissa:.2f}\mathrm{{E}}{{{exponent:+03d}}}$"


def apply_log_scale(ax: Axes) -> None:
    """Log-scale both axes, with minor ticks shorter than major ones."""
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.tick_params(
        which="minor",
        length=plt.rcParams["xtick.major.size"] * MINOR_TICK_LENGTH_RATIO,
    )

    # On these panels matplotlib places minor ticks on one log axis but not the
    # other (measured here: 0 on x, 80 on y). Uneven minors read as a bug, so
    # drop both. The locators answer without a canvas draw.
    x_minors = len(ax.xaxis.get_minor_locator()())
    y_minors = len(ax.yaxis.get_minor_locator()())
    if (x_minors > 0) != (y_minors > 0):
        ax.xaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_minor_locator(NullLocator())


def point_density(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Estimate one 2D density value per point, in the space of the arrays given."""
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        raise ValueError("Cannot estimate point density: x or y has zero range")

    if len(x) <= DENSITY_EXACT_MAX_N:
        values = np.vstack([x, y])
        return gaussian_kde(values)(values)
    return _grid_density(x, y)


def _grid_density(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Fit a subsampled KDE, evaluate it on a grid, interpolate it back to every point."""
    rng = np.random.default_rng(DENSITY_SEED)
    picked = rng.choice(len(x), size=min(len(x), DENSITY_KDE_SUBSAMPLE_N), replace=False)
    kernel = gaussian_kde(np.vstack([x[picked], y[picked]]))

    x_edges = np.linspace(x.min(), x.max(), DENSITY_GRID_SIZE)
    y_edges = np.linspace(y.min(), y.max(), DENSITY_GRID_SIZE)
    x_mesh, y_mesh = np.meshgrid(x_edges, y_edges, indexing="ij")
    grid = kernel(np.vstack([x_mesh.ravel(), y_mesh.ravel()])).reshape(
        DENSITY_GRID_SIZE, DENSITY_GRID_SIZE
    )

    # fill_value=None extrapolates: rounding puts the min/max points a hair
    # outside the grid they defined.
    interpolate = RegularGridInterpolator(
        (x_edges, y_edges), grid, method="linear", bounds_error=False, fill_value=None
    )
    return interpolate(np.column_stack([x, y]))


# =============================================================================
# 6. ASSEMBLY
# =============================================================================
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
def read_panel_tables(config: FigureS1Config) -> tuple[pd.DataFrame, ...]:
    """Read the three parquet tables and fail once if a required column is missing."""
    paths = (config.pbl_pbr, config.tech, config.bio)
    tables = []

    for path, spec in zip(paths, PANEL_SPECS, strict=True):
        table = pd.read_parquet(path)
        missing_columns = {spec.x, spec.y} - set(table.columns)
        if missing_columns:
            raise ValueError(f"{path} is missing required columns: {sorted(missing_columns)}")
        tables.append(table)

    return tuple(tables)


def read_tsv_columns(path: Path, required: tuple[str, ...]) -> pd.DataFrame:
    """Read only the required TSV columns and fail once if any is missing."""
    table = pd.read_csv(path, sep="\t", usecols=lambda column: column in required)
    missing_columns = set(required) - set(table.columns)
    if missing_columns:
        raise ValueError(f"{path} is missing required columns: {sorted(missing_columns)}")
    return table


def run(config: FigureS1Config) -> None:
    """Build Figure S1 end to end: canvas PDF and review raster."""
    config.validate()
    apply_house_style()

    tables = read_panel_tables(config)
    frames = tuple(
        build_scatter_frame(table, spec) for table, spec in zip(tables, PANEL_SPECS, strict=True)
    )
    annotations = read_tsv_columns(config.annotations, ANNOTATION_COLUMNS)
    gene_counts = build_insertion_counts(annotations)
    count_histogram = build_count_histogram(gene_counts)

    layout = build_layout()
    figure = draw_figure(frames, count_histogram, layout)
    # cnsplots sets savefig.bbox to "tight"; pin the page to the designed canvas
    # for the same reason Figure 1 does.
    figure.savefig(
        config.output_pdf,
        format="pdf",
        bbox_inches=Bbox.from_bounds(0.0, 0.0, *figure.get_size_inches()),
    )
    rasterize_review(config.output_pdf, config.output_review)

    logger.success(
        f"Wrote {config.output_pdf} "
        f"({CANVAS_WIDTH_PT:.0f} x {layout.canvas_height_pt:.2f} pt) "
        f"and {config.output_review}"
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Build the composite Figure S1")
    parser.add_argument("--pbl-pbr", type=Path, required=True, help="PBL vs PBR parquet")
    parser.add_argument("--tech", type=Path, required=True, help="Technical replicate parquet")
    parser.add_argument("--bio", type=Path, required=True, help="Biological replicate parquet")
    parser.add_argument(
        "--annotations", type=Path, required=True, help="Panel D insertion annotations TSV"
    )
    parser.add_argument("--output-pdf", type=Path, required=True, help="Final figure PDF")
    parser.add_argument("--output-review", type=Path, required=True, help="Review PNG")
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
    config = FigureS1Config(
        pbl_pbr=args.pbl_pbr,
        tech=args.tech,
        bio=args.bio,
        annotations=args.annotations,
        output_pdf=args.output_pdf,
        output_review=args.output_review,
    )
    try:
        run(config)
    except Exception as e:
        logger.exception(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
