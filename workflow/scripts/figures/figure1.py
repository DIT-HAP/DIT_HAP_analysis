#!/usr/bin/env python3

"""
Figure 1 — DIT-HAP concept, LFC coverage, spike-in linearity
============================================================

Assembles the first main-text figure from three panels:

  A  the DIT-HAP workflow schematic, drawn outside the pipeline and supplied as
     a one-page vector PDF;
  B  LFC gene coverage of Traditional Tn-Seq versus DIT-HAP, by FYPO viability
     class;
  C  spike-in dilution linearity.

Every drawing decision for this figure lives in this file. The script imports
no stage module from ``workflow/src`` and calls no upstream renderer, so
changing Figure 1 cannot change another stage's output.

Panel A is vector for a reason: matplotlib cannot embed PDF or SVG as vector
content (``imshow`` always rasterizes), so the canvas is drawn with A's region
left blank and ``pypdf`` then scales the schematic page into that region and
merges the two pages. The review PNG is rasterized from the merged PDF, not
from the canvas, so it shows the same picture the PDF does.

Input
-----
- results/2a_coverage/cross_dataset/lfc_coverage_comparison.tsv
  Columns: dataset, fypo_class, n_total, n_dr, n_lfc, pct_dr, pct_lfc.
- results/1a_pcr_qc/spike_in_stats.tsv
  Columns: Sample, Name, Relative_Dilution_Ratio, Relative_Read_Ratio.
- resources/figures/figure1a.pdf
  One page; its page box is its content box.

Output
------
- results/7a_figures/figure1.pdf                merged vector figure
- results/7a_figures/figure1.review.png         360 dpi review raster
- results/7a_figures/_work/figure1_canvas.pdf   canvas before the overlay (debug)

Usage
-----
    python workflow/scripts/figures/figure1.py \\
        --lfc-comparison results/2a_coverage/cross_dataset/lfc_coverage_comparison.tsv \\
        --spikein results/1a_pcr_qc/spike_in_stats.tsv \\
        --schematic resources/figures/figure1a.pdf \\
        --output-pdf results/7a_figures/figure1.pdf \\
        --output-review results/7a_figures/figure1.review.png \\
        --canvas results/7a_figures/_work/figure1_canvas.pdf

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-07
Version:  1.0.0
"""

# =============================================================================
# 1. IMPORTS (Standard Library -> Third-Party)
# =============================================================================
import argparse
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
import subprocess
import sys
from typing import cast

import cnsplots as cns
from loguru import logger
from matplotlib.axes import Axes
from matplotlib.container import BarContainer
from matplotlib.figure import Figure
from matplotlib.transforms import Bbox
import pandas as pd
from pypdf import PdfReader, PdfWriter, Transformation

# =============================================================================
# 2. CONSTANTS (house style, layout, data vocabulary)
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
FURNITURE_COLOR = cns.GRAY

# Series colours are palette positions 2/3, the greyscale-separable pair the
# house style reserves for a two-series contrast; panel B's highlight band uses
# position 0, the accent. Resolved from the palette so a palette change
# re-colours the figure.
SERIES_COLORS = cns.get_hexcolors_from_apalette([2, 3], HOUSE_PALETTE)
ACCENT_COLOR = cns.get_hexcolors_from_apalette([0], HOUSE_PALETTE)[0]
FOCUS_BAND_ALPHA = 0.08

# --- Canvas and panel geometry, in points -----------------------------------
# 493 pt is the Cell two-column width. Panel A is the full width, so it sets
# the canvas width; every other box is placed explicitly below.
CANVAS_WIDTH_PT = 493.0
BOTTOM_MARGIN_PT = 2.0
PANEL_LETTER_ROOM_PT = 10.0  # clear band above each row for its letter
PANEL_LETTER_GAP_PT = 2.0
PANEL_LETTER_FONT_PT = 8.0
ROW_GAP_PT = 6.0  # between panel A and the bottom row

SCHEMATIC_WIDTH_PT = CANVAS_WIDTH_PT

# Bottom row. The user tunes these after seeing a render; the values here are
# the first cut from docs/figures.md.
ROW_PANEL_HEIGHT_PT = 80.0
ROW_BELOW_DECORATION_PT = 42.0  # x ticks + xlabel + panel B's legend
PANEL_B_X_PT = 130.0
PANEL_B_WIDTH_PT = 120.0
PANEL_C_X_PT = 300.0
PANEL_C_WIDTH_PT = 110.0

# --- Panel B vocabulary -----------------------------------------------------
COVERAGE_CLASSES = ("inviable", "viable", "depends_on_conditions", "unknown", "total")
CLASS_LABELS = {
    "inviable": "Inviable",
    "viable": "Viable",
    "depends_on_conditions": "Depends on\nconditions",
    "unknown": "Unknown",
    "total": "Total",
}
COVERAGE_DATASETS = ("LD_haploid", "HD_DIT_HAP")
DATASET_LABELS = {"LD_haploid": "Traditional Tn-Seq", "HD_DIT_HAP": "DIT-HAP"}

PERCENT_PER_UNIT = 100.0
AXIS_HEADROOM_PCT = 112
AXIS_TICK_VALUES = [0, 20, 40, 60, 80, 100]
COVERAGE_X_LABEL = "Covered genes (%)"
COVERAGE_Y_LABEL = "FYPO viability class"
FOCUS_CLASS_INDEX = 0  # Inviable is first in COVERAGE_CLASSES
BAR_COUNT_FONT_PT = 5  # fit-critical: the count must clear an 80 pt panel
LEGEND_BELOW_Y = -0.34

# --- Panel C vocabulary -----------------------------------------------------
SPIKEIN_REFERENCE_SAMPLE = "Spikein0"
SPIKEIN_SITE_PREFIX = "Spike-in Insertion "
SPIKEIN_SITE_LABEL_PREFIX = "#"
SPIKEIN_X_LABEL = r"$\log_2$(relative dilution ratio)"
SPIKEIN_Y_LABEL = r"$\log_2$(relative read ratio)"

# --- Input contracts --------------------------------------------------------
COVERAGE_COLUMNS = ("dataset", "fypo_class", "n_lfc", "pct_lfc")
SPIKEIN_COLUMNS = ("Sample", "Name", "Relative_Dilution_Ratio", "Relative_Read_Ratio")


# =============================================================================
# 3. CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class Figure1Config:
    """Resolved input/output paths for Figure 1."""

    lfc_comparison: Path
    spikein: Path
    schematic: Path
    output_pdf: Path
    output_review: Path
    canvas_pdf: Path

    def validate(self) -> None:
        """Raise if a required input is missing, then create the output dirs."""
        for path in [self.lfc_comparison, self.spikein, self.schematic]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for output in [self.output_pdf, self.output_review, self.canvas_pdf]:
            output.parent.mkdir(parents=True, exist_ok=True)


@dataclass(kw_only=True)
class Layout:
    """Resolved point geometry, derived from the schematic's page aspect."""

    canvas_height_pt: float
    schematic_height_pt: float
    row_top_pt: float
    schematic_rect: tuple[float, float, float, float]
    panel_b_rect: tuple[float, float, float, float]
    panel_c_rect: tuple[float, float, float, float]


# =============================================================================
# 4. DATA PREPARATION (pure; schemas are validated at the I/O boundary)
# =============================================================================
def build_coverage_frame(comparison: pd.DataFrame) -> pd.DataFrame:
    """Slice the two libraries into the long bar form, in display order."""
    missing_classes = set(COVERAGE_CLASSES) - set(comparison["fypo_class"])
    if missing_classes:
        raise ValueError(f"LFC comparison is missing FYPO classes: {sorted(missing_classes)}")

    is_selected = comparison["dataset"].isin(COVERAGE_DATASETS)
    frame = comparison.loc[is_selected, ["dataset", "fypo_class", "n_lfc", "pct_lfc"]].copy()
    frame["coverage_pct"] = frame["pct_lfc"] * PERCENT_PER_UNIT
    frame["class_label"] = frame["fypo_class"].map(CLASS_LABELS)
    frame["dataset_label"] = frame["dataset"].map(DATASET_LABELS)
    return frame


def build_spikein_frame(spikein: pd.DataFrame) -> pd.DataFrame:
    """Drop the zero-dilution reference and label each site as #1..#5."""
    is_measured = spikein["Sample"] != SPIKEIN_REFERENCE_SAMPLE
    frame = spikein.loc[is_measured].copy()
    frame["Name"] = frame["Name"].str.replace(SPIKEIN_SITE_PREFIX, SPIKEIN_SITE_LABEL_PREFIX)
    return frame.sort_values("Name", kind="stable")


def build_layout(schematic_width_pt: float, schematic_height_pt: float) -> Layout:
    """Derive the canvas height and every axes box from the schematic's aspect."""
    schematic_height = SCHEMATIC_WIDTH_PT * schematic_height_pt / schematic_width_pt
    row_top = PANEL_LETTER_ROOM_PT + schematic_height + ROW_GAP_PT + PANEL_LETTER_ROOM_PT
    canvas_height = row_top + ROW_PANEL_HEIGHT_PT + ROW_BELOW_DECORATION_PT + BOTTOM_MARGIN_PT

    return Layout(
        canvas_height_pt=canvas_height,
        schematic_height_pt=schematic_height,
        row_top_pt=row_top,
        schematic_rect=_axes_rect(
            0.0, PANEL_LETTER_ROOM_PT, SCHEMATIC_WIDTH_PT, schematic_height, canvas_height
        ),
        panel_b_rect=_axes_rect(
            PANEL_B_X_PT, row_top, PANEL_B_WIDTH_PT, ROW_PANEL_HEIGHT_PT, canvas_height
        ),
        panel_c_rect=_axes_rect(
            PANEL_C_X_PT, row_top, PANEL_C_WIDTH_PT, ROW_PANEL_HEIGHT_PT, canvas_height
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


def draw_figure(coverage: pd.DataFrame, spikein: pd.DataFrame, layout: Layout) -> Figure:
    """Draw the canvas: panels B and C, with panel A's region left blank."""
    fig = cns.figure(width=CANVAS_WIDTH_PT, height=layout.canvas_height_pt, unit="pt")
    schematic_ax = fig.add_axes(layout.schematic_rect)
    schematic_ax.set_axis_off()

    coverage_ax = fig.add_axes(layout.panel_b_rect)
    render_coverage_panel(coverage_ax, coverage)
    spikein_ax = fig.add_axes(layout.panel_c_rect)
    render_spikein_panel(spikein_ax, spikein)

    place_panel_letters(
        fig,
        (
            (schematic_ax, "A", PANEL_LETTER_ROOM_PT),
            (coverage_ax, "B", layout.row_top_pt),
            (spikein_ax, "C", layout.row_top_pt),
        ),
        layout.canvas_height_pt,
    )
    return fig


def place_panel_letters(
    fig: Figure,
    letters: Sequence[tuple[Axes, str, float]],
    canvas_height_pt: float,
) -> None:
    """Place each panel letter just above and left of its panel's y-axis decoration.

    The x comes from the axes' tight bbox, so the letter clears the y tick labels
    and the ylabel whatever they measure rather than sitting on the spine.
    """
    fig.canvas.draw()

    for ax, label, y_top_pt in letters:
        bbox = ax.get_tightbbox()
        assert bbox is not None
        left_pt = max(0.0, bbox.x0 * 72 / fig.dpi - PANEL_LETTER_GAP_PT)
        fig.text(
            left_pt / CANVAS_WIDTH_PT,
            1.0 - (y_top_pt - PANEL_LETTER_GAP_PT) / canvas_height_pt,
            label,
            ha="left",
            va="bottom",
            fontweight="bold",
            fontsize=PANEL_LETTER_FONT_PT,
        )


def render_coverage_panel(ax: Axes, frame: pd.DataFrame) -> None:
    """Draw grouped horizontal bars: coverage percent per FYPO class, two libraries."""
    class_labels = [CLASS_LABELS[cls] for cls in COVERAGE_CLASSES]
    cns.barplot(
        data=frame,
        x="coverage_pct",
        y="class_label",
        orient="h",
        order=class_labels,
        hue="dataset_label",
        hue_order=[DATASET_LABELS[dataset] for dataset in COVERAGE_DATASETS],
        palette=SERIES_COLORS,
        ax=ax,
    )

    ax.axhspan(
        FOCUS_CLASS_INDEX - 0.5,
        FOCUS_CLASS_INDEX + 0.5,
        color=ACCENT_COLOR,
        alpha=FOCUS_BAND_ALPHA,
        linewidth=0,
        zorder=0,
    )
    focus_tick = ax.get_yticklabels()[FOCUS_CLASS_INDEX]
    focus_tick.set_color(ACCENT_COLOR)
    focus_tick.set_fontweight("bold")

    add_count_labels(ax, frame)
    ax.set_xlim(0, AXIS_HEADROOM_PCT)
    ax.set_xticks(AXIS_TICK_VALUES)
    ax.set_xlabel(COVERAGE_X_LABEL)
    ax.set_ylabel(COVERAGE_Y_LABEL)

    handles, legend_labels = ax.get_legend_handles_labels()
    ax.legend(
        handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, LEGEND_BELOW_Y),
        ncol=len(COVERAGE_DATASETS),
        frameon=False,
    )


def add_count_labels(ax: Axes, frame: pd.DataFrame) -> None:
    """Label every bar with its absolute covered-gene count."""
    for dataset, container in zip(COVERAGE_DATASETS, ax.containers, strict=True):
        counts = (
            frame[frame["dataset"] == dataset]
            .set_index("fypo_class")
            .loc[list(COVERAGE_CLASSES), "n_lfc"]
        )
        labels = [f"{int(count):,}" for count in counts]
        ax.bar_label(
            cast(BarContainer, container),
            labels=labels,
            padding=1,
            fontsize=BAR_COUNT_FONT_PT,
            color=FURNITURE_COLOR,
        )


def render_spikein_panel(ax: Axes, frame: pd.DataFrame) -> None:
    """Draw the dilution-linearity scatter with one overall fit."""
    # color=<column> (not hue=) gives one legend entry per site but a single
    # overall fit; hue= would fit and annotate five separate lines.
    cns.regplot(
        data=frame,
        x="Relative_Dilution_Ratio",
        y="Relative_Read_Ratio",
        color="Name",
        add_equation=True,
        ax=ax,
    )
    ax.set_xlabel(SPIKEIN_X_LABEL)
    ax.set_ylabel(SPIKEIN_Y_LABEL)

    # cnsplots parks the colour column's name as the legend title; replace it.
    cns.take_legend_out(title="Spike-in", ax=ax)


# =============================================================================
# 6. ASSEMBLY (vector overlay and review raster)
# =============================================================================
def page_size_pt(pdf: Path) -> tuple[float, float]:
    """Return the first page's (width, height) in points."""
    page = PdfReader(str(pdf)).pages[0]
    return float(page.mediabox.width), float(page.mediabox.height)


def overlay_schematic(
    canvas_pdf: Path, schematic_pdf: Path, output_pdf: Path, layout: Layout
) -> None:
    """Scale the schematic uniformly into panel A's box and merge it onto the canvas."""
    writer = PdfWriter(clone_from=str(canvas_pdf))
    canvas_page = writer.pages[0]
    schematic_page = PdfReader(str(schematic_pdf)).pages[0]

    source_width = float(schematic_page.mediabox.width)
    source_height = float(schematic_page.mediabox.height)
    scale = min(SCHEMATIC_WIDTH_PT / source_width, layout.schematic_height_pt / source_height)

    x_pt = (SCHEMATIC_WIDTH_PT - source_width * scale) / 2
    y_top_to_bottom = layout.canvas_height_pt - PANEL_LETTER_ROOM_PT - layout.schematic_height_pt
    y_pt = y_top_to_bottom + (layout.schematic_height_pt - source_height * scale) / 2

    canvas_page.merge_transformed_page(
        schematic_page, Transformation().scale(scale).translate(x_pt, y_pt)
    )
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
def read_tsv(path: Path, required_columns: tuple[str, ...]) -> pd.DataFrame:
    """Read a TSV and fail once if any required column is missing."""
    table = pd.read_csv(path, sep="\t")
    missing_columns = set(required_columns) - set(table.columns)
    if missing_columns:
        raise ValueError(f"{path} is missing required columns: {sorted(missing_columns)}")
    return table


def run(config: Figure1Config) -> None:
    """Build Figure 1 end to end: canvas, vector overlay, review raster."""
    config.validate()
    apply_house_style()

    coverage = build_coverage_frame(read_tsv(config.lfc_comparison, COVERAGE_COLUMNS))
    spikein = build_spikein_frame(read_tsv(config.spikein, SPIKEIN_COLUMNS))

    schematic_width, schematic_height = page_size_pt(config.schematic)
    layout = build_layout(schematic_width, schematic_height)

    figure = draw_figure(coverage, spikein, layout)
    # cnsplots sets savefig.bbox to "tight", which would crop the page to the
    # drawn artists and break the point-exact overlay below. Save the designed
    # canvas instead.
    figure.savefig(
        config.canvas_pdf,
        format="pdf",
        bbox_inches=Bbox.from_bounds(0.0, 0.0, *figure.get_size_inches()),
    )
    overlay_schematic(config.canvas_pdf, config.schematic, config.output_pdf, layout)
    rasterize_review(config.output_pdf, config.output_review)

    logger.success(
        f"Wrote {config.output_pdf} ({CANVAS_WIDTH_PT:.0f} x {layout.canvas_height_pt:.2f} pt) "
        f"and {config.output_review}"
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Build the composite Figure 1")
    parser.add_argument(
        "--lfc-comparison", type=Path, required=True, help="Combined LFC comparison TSV"
    )
    parser.add_argument("--spikein", type=Path, required=True, help="Spike-in stats TSV")
    parser.add_argument("--schematic", type=Path, required=True, help="Panel A schematic PDF")
    parser.add_argument("--output-pdf", type=Path, required=True, help="Final figure PDF")
    parser.add_argument("--output-review", type=Path, required=True, help="Review PNG")
    parser.add_argument("--canvas", type=Path, required=True, help="Canvas PDF before the overlay")
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
    config = Figure1Config(
        lfc_comparison=args.lfc_comparison,
        spikein=args.spikein,
        schematic=args.schematic,
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
