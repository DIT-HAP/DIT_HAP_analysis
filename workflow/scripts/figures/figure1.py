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

The page is a ``cns.multipanel`` flow: A's reserved box fills the first row,
``mp.newline()`` puts B and C at the left of the second, and the unused right
half of that row stays blank. Each panel's box is its axes size plus whatever
y-axis decoration multipanel measures at render time, so B and C keep their own
labels clear instead of overlapping; the letters are placed by multipanel at the
left edge of each panel's column.

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
# 493 pt is the Cell two-column width, which multipanel takes as the page width.
CANVAS_WIDTH_PT = 493.0
ROW_GAP_PT = 6.0  # between panel A and the bottom row

# Panel A alone spans the page. multipanel reserves each panel's letter column
# on the left, and A is the only panel that must reach the right edge, so its
# axes width gives that column back: 6 pt is the rendered width of a bold 8 pt
# "A". The schematic is then scaled into the axes box after the page is saved.
SCHEMATIC_LABEL_COLUMN_PT = 6.0
SCHEMATIC_WIDTH_PT = CANVAS_WIDTH_PT - SCHEMATIC_LABEL_COLUMN_PT
# multipanel's pad_top is in display pixels (figure_dpi / 72 per point), so the
# letter band above A is the label's own height plus this clearance.
PANEL_LETTER_PAD_TOP_PX = round(2.8 * HOUSE_FIGURE_DPI / 72)

# Bottom row.
ROW_PANEL_HEIGHT_PT = 80.0
ROW_BELOW_DECORATION_PT = 42.0  # x ticks + xlabel + panel B's legend
PANEL_B_WIDTH_PT = 120.0
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
    coverage: pd.DataFrame, spikein: pd.DataFrame, schematic_height_pt: float
) -> tuple[cns.multipanel, Axes]:
    """Draw the page: A's reserved schematic box, then B and C left-aligned below.

    Returns the multipanel itself, not just its figure: cnsplots keeps the
    measured-layout draw callback as a weak reference, so the caller must hold
    the object until the page is written or the panels stay unmeasured.
    """
    mp = cns.multipanel(max_width=CANVAS_WIDTH_PT, unit="pt")

    schematic_ax = mp.panel(
        "A",
        width=SCHEMATIC_WIDTH_PT,
        height=schematic_height_pt,
        pad_top=PANEL_LETTER_PAD_TOP_PX,
        margin_bottom=ROW_GAP_PT,
    )
    hide_schematic_axes(schematic_ax)

    mp.newline()  # starts B's row at the page's left edge, under A
    coverage_ax = mp.panel(
        "B",
        width=PANEL_B_WIDTH_PT,
        height=ROW_PANEL_HEIGHT_PT,
        margin_bottom=ROW_BELOW_DECORATION_PT,
    )
    render_coverage_panel(coverage_ax, coverage)
    spikein_ax = mp.panel(
        "C",
        width=PANEL_C_WIDTH_PT,
        height=ROW_PANEL_HEIGHT_PT,
        margin_bottom=ROW_BELOW_DECORATION_PT,
    )
    render_spikein_panel(spikein_ax, spikein)
    return mp, schematic_ax


def hide_schematic_axes(ax: Axes) -> None:
    """Blank panel A: its content arrives as vector PDF after the page is saved.

    The empty tick sets matter: multipanel measures y-axis decoration even on an
    axis-off axes, and the placeholder ticks would widen A's left reserve until
    its box ran off the right edge of the page.
    """
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()


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


def axes_box_pt(ax: Axes, figure: Figure) -> tuple[float, float, float, float]:
    """Return an axes box as (x, y from the page top, width, height) in points."""
    width_in, height_in = figure.get_size_inches()
    width_pt, height_pt = width_in * 72, height_in * 72
    box = ax.get_position()
    return box.x0 * width_pt, (1 - box.y1) * height_pt, box.width * width_pt, box.height * height_pt


def overlay_schematic(
    canvas_pdf: Path,
    schematic_pdf: Path,
    output_pdf: Path,
    box_pt: tuple[float, float, float, float],
) -> None:
    """Scale the schematic uniformly into panel A's box and merge it onto the canvas."""
    x_pt, y_top_pt, box_width_pt, box_height_pt = box_pt
    _, canvas_height_pt = page_size_pt(canvas_pdf)
    writer = PdfWriter(clone_from=str(canvas_pdf))
    canvas_page = writer.pages[0]
    schematic_page = PdfReader(str(schematic_pdf)).pages[0]

    source_width = float(schematic_page.mediabox.width)
    source_height = float(schematic_page.mediabox.height)
    scale = min(box_width_pt / source_width, box_height_pt / source_height)

    x_pt += (box_width_pt - source_width * scale) / 2
    y_from_top = y_top_pt + (box_height_pt - source_height * scale) / 2
    # PDF coordinates count from the page's bottom edge.
    y_pt = canvas_height_pt - y_from_top - source_height * scale

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


def run(config: Figure1Config) -> tuple[float, float]:
    """Build Figure 1 end to end and return the page's (width, height) in points."""
    config.validate()
    apply_house_style()

    coverage = build_coverage_frame(read_tsv(config.lfc_comparison, COVERAGE_COLUMNS))
    spikein = build_spikein_frame(read_tsv(config.spikein, SPIKEIN_COLUMNS))

    schematic_width, schematic_height = page_size_pt(config.schematic)
    schematic_box_height = SCHEMATIC_WIDTH_PT * schematic_height / schematic_width

    mp, schematic_ax = draw_figure(coverage, spikein, schematic_box_height)
    assert mp.fig is not None
    # Settle multipanel's measured layout, then save the designed canvas:
    # cnsplots sets savefig.bbox to "tight", which would crop the page to the
    # drawn artists and shift panel A's box away from the overlay below.
    mp.fig.canvas.draw()
    cns.savefig(config.canvas_pdf, fig=mp.fig, bbox_inches=None)
    overlay_schematic(
        config.canvas_pdf,
        config.schematic,
        config.output_pdf,
        axes_box_pt(schematic_ax, mp.fig),
    )
    rasterize_review(config.output_pdf, config.output_review)

    width_pt, height_pt = mp.fig.get_size_inches() * 72
    logger.success(
        f"Wrote {config.output_pdf} ({width_pt:.0f} x {height_pt:.2f} pt) "
        f"and {config.output_review}"
    )
    return width_pt, height_pt


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
