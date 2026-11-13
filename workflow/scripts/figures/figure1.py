#!/usr/bin/env python3

"""
Figure 1 — DIT-HAP concept, coverage, spike-in, LFC trajectories, FACS
======================================================================

Assembles the first main-text figure from five panels:

  A  the DIT-HAP workflow schematic, drawn outside the pipeline and supplied as
     a one-page vector PDF;
  B  LFC gene coverage of Traditional Tn-Seq versus DIT-HAP, by FYPO viability
     class;
  C  spike-in dilution linearity;
  D  gene-level LFC through the 12.5-generation time course, one trajectory
     per gene, coloured by the two FYPO viability classes B highlights;
  E  FACS-validated genes: each gene's FACS depletion curves (two independent
     strains) against its DIT-HAP gene-level LFC, one small axes per gene.

Every drawing decision for this figure lives in this file. The script imports
no stage module from ``workflow/src`` and calls no upstream renderer, so
changing Figure 1 cannot change another stage's output.

Panel A is vector for a reason: matplotlib cannot embed PDF or SVG as vector
content (``imshow`` always rasterizes), so the canvas is drawn with A's region
left blank and ``pypdf`` then copies the schematic page into that region at
scale 1.0 (the panel box is the schematic's own page box) and merges the two
pages. The review PNG is rasterized from the merged PDF, not from the canvas,
so it shows the same picture the PDF does.

The page is an A4 manuscript page (210 x 297 mm). ``cns.multipanel`` sizes only
the canvas width and lets the rows decide the height, so the flow below draws a
content-sized canvas and the assembly mounts it 40 pt in from the A4 page's
top-left corner, leaving the rest of the page blank. Within the content, A takes
the schematic page's actual size at the left of the first row,
the first ``mp.newline()`` puts B, C and D at the left of the second, and the
second puts E across the third. Each panel's box is its axes size plus whatever
y-axis decoration multipanel measures at render time, so a wide ylabel pushes
only its own panel right; the letters are placed by multipanel at the left edge
of each panel's column.

Input
-----
- results/2a_coverage/cross_dataset/lfc_coverage_comparison.tsv
  Columns: dataset, fypo_class, n_total, n_dr, n_lfc, pct_dr, pct_lfc.
- results/1a_pcr_qc/spike_in_stats.tsv
  Columns: Sample, Name, Relative_Dilution_Ratio, Relative_Read_Ratio.
- <HD_DIT_HAP release>/gene_level/LFC.tsv
  Columns: Systematic ID, Name, FYPOviability, YES0..YES4 (log2 fold change
  per generation; YES0 is the 0-generation baseline).
- resources/curated/FACS_depletion_curves_from_XH.xlsx
  Columns: Name, Strain_alias, Generations, LFC (one row per gene x strain x
  timepoint; FACS LFC is depletion-positive, the opposite sign to the release).
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
        --gene-lfc <HD_DIT_HAP release>/gene_level/LFC.tsv \\
        --facs resources/curated/FACS_depletion_curves_from_XH.xlsx \\
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
from math import ceil, floor
from pathlib import Path
import subprocess
import sys
from typing import cast

import cnsplots as cns
from loguru import logger
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.container import BarContainer
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
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
# The manuscript page is A4 (210 x 297 mm). multipanel sizes the figure by
# width only, so the content canvas is A4 wide minus the page margins and the
# assembly mounts it PAGE_MARGIN_PT in from the page's top-left corner; the
# rest of the page stays blank.
MM_PER_INCH = 25.4
A4_WIDTH_PT = 210.0 / MM_PER_INCH * 72.0  # 595.28
A4_HEIGHT_PT = 297.0 / MM_PER_INCH * 72.0  # 841.89
PAGE_MARGIN_PT = 40.0  # white border at the left, top and right
# The content width is the page width minus the left and right margins, so the
# right margin comes for free and the figure cannot run off the page.
CANVAS_WIDTH_PT = A4_WIDTH_PT - 2 * PAGE_MARGIN_PT  # 515.28
ROW_GAP_PT = 6.0  # between panel rows

# Panel E spans the content width. Panel A instead takes the schematic PDF's own
# page box, so the artwork is embedded at its actual size (scale 1.0) rather
# than blown up to the page width. multipanel reserves each panel's letter
# column on the left, so a page-spanning axis gives that column back: 6 pt is
# the rendered width of a bold 8 pt "A".
PANEL_LABEL_COLUMN_PT = 6.0
PAGE_CONTENT_WIDTH_PT = CANVAS_WIDTH_PT - PANEL_LABEL_COLUMN_PT
# multipanel's pad_top is in display pixels (figure_dpi / 72 per point), so the
# letter band above A is the label's own height plus this clearance.
PANEL_LETTER_PAD_TOP_PX = round(2.8 * HOUSE_FIGURE_DPI / 72)

# Bottom row: B, C and D share one square footprint so the three plotting areas
# read as a set; 100 pt is the house SQUARE baseline (inside the 80-150 pt band)
# and the three columns still fit the page width with their own y decorations.
ROW_PANEL_SIZE_PT = 100.0
# x ticks + xlabel; the legends sit above the row, and multipanel measures
# them as top decoration.
ROW_BELOW_DECORATION_PT = 26.0

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
COVERAGE_LEGEND_TITLE = "Library"
DATASET_LABELS = {"LD_haploid": "Traditional Tn-Seq", "HD_DIT_HAP": "DIT-HAP"}

PERCENT_PER_UNIT = 100.0
AXIS_HEADROOM_PCT = 112
AXIS_TICK_VALUES = [0, 20, 40, 60, 80, 100]
COVERAGE_X_LABEL = "Covered genes (%)"
COVERAGE_Y_LABEL = "FYPO viability class"
FOCUS_CLASS_INDEX = 0  # Inviable is first in COVERAGE_CLASSES
BAR_COUNT_FONT_PT = 5  # fit-critical: the count must clear inside its bar
# A data hairline would vanish as a legend swatch, so the legend draws its own.
LEGEND_SWATCH_LINE_WIDTH_PT = 0.8
# Narrower than matplotlib's 2 em default: the horizontal legends are wide and
# sit over their own panel only.
LEGEND_COLUMN_SPACING_EM = 0.8

# --- Panel C vocabulary -----------------------------------------------------
SPIKEIN_REFERENCE_SAMPLE = "Spikein0"
SPIKEIN_SITE_PREFIX = "Spike-in Insertion "
SPIKEIN_SITE_LABEL_PREFIX = "#"
SPIKEIN_LEGEND_TITLE = "Spike-in"
SPIKEIN_POINT_SIZE = 25.0
SPIKEIN_X_LABEL = r"$\log_2$(relative dilution ratio)"
SPIKEIN_Y_LABEL = r"$\log_2$(relative read ratio)"

# --- Panel D vocabulary -----------------------------------------------------
# The release's five sampling points: the LFC columns in order, and the
# generations they stand for (upstream config's time_points).
GENERATION_COLUMNS = ("YES0", "YES1", "YES2", "YES3", "YES4")
GENERATIONS = (0.0, 2.352, 5.588, 9.104, 12.48)
# Display rounding only (the full-precision values are in the data and docs):
# three decimals collide under a 100 pt axes.
GENERATION_TICK_LABELS = ("0", "2.4", "5.6", "9.1", "12.5")
TRAJECTORY_CLASSES = ("inviable", "viable")
TRAJECTORY_LABELS = {"inviable": "Inviable", "viable": "Viable"}
# Inviable keeps panel B's accent, the colour B already uses for that class;
# viable takes palette slot 1, which no other series in this figure uses.
TRAJECTORY_COLORS = {
    "inviable": ACCENT_COLOR,
    "viable": cns.get_hexcolors_from_apalette([1], HOUSE_PALETTE)[0],
}
# Thousands of gene curves share one 100 pt axes, so each is a hairline and
# transparent: dense regions darken by overplotting, single genes stay legible.
TRAJECTORY_LINE_WIDTH_PT = 0.5
TRAJECTORY_LINE_ALPHA = 0.05
# Inviable draws over viable so the highlighted class keeps its accent where
# the two overlap (at the shared 0-generation baseline).
TRAJECTORY_ZORDER = {"inviable": 3, "viable": 2}
LFC_X_LABEL = "Generations"
LFC_Y_LABEL = r"Gene-level LFC ($\log_2$)"
# D's legend straddles its axes' top edge, so only its upper half counts as
# top decoration; D's pad_top gives the rest back to keep D on B/C's baseline.
LFC_LEGEND_HALF_HEIGHT_PT = 8.4

# --- Panel E vocabulary -----------------------------------------------------
DIT_HAP_SERIES = "DIT-HAP"
FACS_SERIES_PREFIX = "FACS "
# The workbook covers ten genes; orc6 and srb4 were retired from the panel by
# user decision (2026-11-12), leaving a clean 4 x 2 grid.
FACS_EXCLUDED_GENES = ("orc6", "srb4")
# The two FACS strains are independent isolates, not replicates, so both are
# drawn and each gets its own Cell-palette colour.
FACS_STRAIN_ORDER = ("Strain#1", "Strain#2")
FACS_STRAIN_COLORS = dict(
    zip(
        FACS_STRAIN_ORDER,
        cns.get_hexcolors_from_apalette([5, 6], HOUSE_PALETTE),
        strict=True,
    )
)
# FACS LFC is depletion-positive (the reference outgrows the mutant); the
# release reports depletion-negative values. Flipping FACS puts both curves on
# panel D's convention, so a depleted gene descends in both.
FACS_LFC_SIGN = -1.0
FACS_X_LABEL = "Generations"
FACS_Y_LABEL = r"LFC ($\log_2$)"
FACS_X_TICKS = (0.0, 4.0, 8.0, 12.0)
FACS_X_LIMITS = (-0.3, 12.8)
FACS_LINE_WIDTH_PT = 2.0  # user-specified, heavier than the 1.5 house default
# A light alpha keeps the three curves from reading as opaque tubes where they
# overlap; the legend swatches stay fully opaque.
FACS_LINE_ALPHA = 0.8
FACS_Y_LABEL_CENTER_X_PT = 15.0  # centre of the rotated shared y label
FACS_X_LABEL_BOTTOM_PT = 3.0  # bottom of the shared x label

# --- Panel E geometry -------------------------------------------------------
# Eight genes in a 4 x 2 grid; the legend is one row above the panel. The mini
# axes are insets of one page-wide container panel, so they follow its box at
# render time.
FACS_GRID_COLUMNS = 4
FACS_MINI_AXES_WIDTH_PT = 80.0
FACS_MINI_AXES_HEIGHT_PT = 60.0
FACS_COLUMN_GAP_PT = 10.0
FACS_ROW_GAP_PT = 16.0  # clears the lower row's 8 pt gene title
FACS_TITLE_BAND_PT = 13.0  # above the upper row's gene titles
FACS_BOTTOM_DECORATION_PT = 22.0  # x tick labels + the shared x label
FACS_LEFT_DECORATION_PT = 34.0  # y tick labels + the shared y label
# DIT-HAP takes the green of slot 4; the two FACS strains take slots 5 and 6
# (salmon, purple). No other panel uses those three slots.
DIT_HAP_COLOR = cns.get_hexcolors_from_apalette([4], HOUSE_PALETTE)[0]

# --- Input contracts --------------------------------------------------------
COVERAGE_COLUMNS = ("dataset", "fypo_class", "n_lfc", "pct_lfc")
SPIKEIN_COLUMNS = ("Sample", "Name", "Relative_Dilution_Ratio", "Relative_Read_Ratio")
TRAJECTORY_COLUMNS = ("Systematic ID", "FYPOviability", *GENERATION_COLUMNS)
GENE_LFC_COLUMNS = ("Systematic ID", "Name", "FYPOviability", *GENERATION_COLUMNS)
FACS_COLUMNS = ("Name", "Strain_alias", "Generations", "LFC")


# =============================================================================
# 3. CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class Figure1Config:
    """Resolved input/output paths for Figure 1."""

    lfc_comparison: Path
    spikein: Path
    gene_lfc: Path
    facs: Path
    schematic: Path
    output_pdf: Path
    output_review: Path
    canvas_pdf: Path

    def validate(self) -> None:
        """Raise if a required input is missing, then create the output dirs."""
        inputs = [self.lfc_comparison, self.spikein, self.gene_lfc, self.facs, self.schematic]
        for path in inputs:
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


def build_trajectory_frame(gene_lfc: pd.DataFrame) -> pd.DataFrame:
    """Reshape the wide per-generation LFC table into one row per gene and generation."""
    is_traced = gene_lfc["FYPOviability"].isin(TRAJECTORY_CLASSES)
    wide = gene_lfc.loc[is_traced, list(TRAJECTORY_COLUMNS)]
    long = wide.melt(
        id_vars=["Systematic ID", "FYPOviability"],
        var_name="generation_column",
        value_name="lfc",
    )
    long["generation"] = long["generation_column"].map(
        dict(zip(GENERATION_COLUMNS, GENERATIONS, strict=True))
    )
    return long


def build_facs_comparison_frame(facs: pd.DataFrame, gene_lfc: pd.DataFrame) -> pd.DataFrame:
    """Stack the FACS and DIT-HAP curves per gene on one sign convention."""
    is_shown = ~facs["Name"].isin(FACS_EXCLUDED_GENES)
    shown_facs = facs.loc[is_shown, list(FACS_COLUMNS)]
    genes = sorted(shown_facs["Name"].unique())
    missing_genes = set(genes) - set(gene_lfc["Name"])
    if missing_genes:
        raise ValueError(f"DIT-HAP LFC table is missing FACS genes: {sorted(missing_genes)}")

    facs_long = shown_facs.rename(columns={"Name": "gene", "Generations": "generation"})
    facs_long["series"] = FACS_SERIES_PREFIX + facs_long["Strain_alias"]
    facs_long["lfc"] = FACS_LFC_SIGN * facs_long["LFC"]

    is_facs_gene = gene_lfc["Name"].isin(genes)
    dit_hap = gene_lfc.loc[is_facs_gene, ["Name", *GENERATION_COLUMNS]]
    dit_hap_long = dit_hap.melt(id_vars="Name", var_name="generation_column", value_name="lfc")
    dit_hap_long["generation"] = dit_hap_long["generation_column"].map(
        dict(zip(GENERATION_COLUMNS, GENERATIONS, strict=True))
    )
    dit_hap_long = dit_hap_long.rename(columns={"Name": "gene"})
    dit_hap_long["series"] = DIT_HAP_SERIES

    columns = ["gene", "series", "generation", "lfc"]
    stacked = pd.concat([dit_hap_long[columns], facs_long[columns]], ignore_index=True)
    return stacked.sort_values(["gene", "series", "generation"], kind="stable")


def facs_gene_rows(frame: pd.DataFrame) -> int:
    """Return the number of grid rows the FACS gene set needs."""
    return ceil(frame["gene"].nunique() / FACS_GRID_COLUMNS)


def facs_panel_height_pt(n_rows: int) -> float:
    """Return panel E's container height for a grid of n_rows rows."""
    axes_height = n_rows * FACS_MINI_AXES_HEIGHT_PT
    gaps_height = (n_rows - 1) * FACS_ROW_GAP_PT
    return FACS_TITLE_BAND_PT + axes_height + gaps_height + FACS_BOTTOM_DECORATION_PT


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
    coverage: pd.DataFrame,
    spikein: pd.DataFrame,
    trajectory: pd.DataFrame,
    facs_comparison: pd.DataFrame,
    schematic_size_pt: tuple[float, float],
) -> tuple[cns.multipanel, Axes]:
    """Draw the page: A's actual-size schematic box, B, C and D, then E below them.

    Returns the multipanel itself, not just its figure: cnsplots keeps the
    measured-layout draw callback as a weak reference, so the caller must hold
    the object until the page is written or the panels stay unmeasured.
    """
    mp = cns.multipanel(max_width=CANVAS_WIDTH_PT, unit="pt")

    schematic_width_pt, schematic_height_pt = schematic_size_pt
    schematic_ax = mp.panel(
        "A",
        width=schematic_width_pt,
        height=schematic_height_pt,
        pad_top=PANEL_LETTER_PAD_TOP_PX,
        margin_bottom=ROW_GAP_PT,
    )
    blank_panel_axes(schematic_ax)

    mp.newline()  # starts B's row at the page's left edge, under A
    coverage_ax = mp.panel(
        "B",
        width=ROW_PANEL_SIZE_PT,
        height=ROW_PANEL_SIZE_PT,
        margin_bottom=ROW_BELOW_DECORATION_PT,
    )
    render_coverage_panel(coverage_ax, coverage)
    spikein_ax = mp.panel(
        "C",
        width=ROW_PANEL_SIZE_PT,
        height=ROW_PANEL_SIZE_PT,
        margin_bottom=ROW_BELOW_DECORATION_PT,
    )
    render_spikein_panel(spikein_ax, spikein)
    lfc_ax = mp.panel(
        "D",
        width=ROW_PANEL_SIZE_PT,
        height=ROW_PANEL_SIZE_PT,
        margin_bottom=ROW_BELOW_DECORATION_PT,
        pad_top=round(LFC_LEGEND_HALF_HEIGHT_PT * HOUSE_FIGURE_DPI / 72),
    )
    render_lfc_trajectory_panel(lfc_ax, trajectory)

    mp.newline()  # panel E spans the page under the B/C/D row
    facs_ax = mp.panel(
        "E",
        width=PAGE_CONTENT_WIDTH_PT,
        height=facs_panel_height_pt(facs_gene_rows(facs_comparison)),
        margin_top=ROW_GAP_PT,
        margin_bottom=0,  # last row: its x label is the page's bottom decoration
    )
    blank_panel_axes(facs_ax)
    render_facs_panel(facs_ax, facs_comparison)
    return mp, schematic_ax


def blank_panel_axes(ax: Axes) -> None:
    """Blank a container panel whose content is insets or a later vector overlay.

    The empty y tick set is what matters: multipanel measures y-axis decoration
    even on an axis-off axes, and the placeholder ticks would widen the panel's
    left reserve until its box ran off the right edge of the page.
    """
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
    add_unframed_legend(ax, columns=len(COVERAGE_DATASETS), title=COVERAGE_LEGEND_TITLE)


def legend_entries(
    ax: Axes, handles: Sequence[Artist] | None = None
) -> tuple[Sequence[Artist], list[str]]:
    """Return a legend's handles and labels, from explicit handles or the axes."""
    if handles is None:
        return ax.get_legend_handles_labels()
    return handles, [cast(str, handle.get_label()) for handle in handles]


def add_unframed_legend(
    ax: Axes,
    columns: int,
    title: str | None = None,
    handles: Sequence[Artist] | None = None,
    loc: str = "lower left",
    anchor: tuple[float, float] = (0.0, 1.0),
) -> None:
    """Draw an unframed legend at an axes-fraction anchor, above the axes by default."""
    entries, labels = legend_entries(ax, handles)
    ax.legend(
        entries,
        labels,
        title=title,
        loc=loc,
        bbox_to_anchor=anchor,
        borderaxespad=0,
        borderpad=0,
        ncol=columns,
        columnspacing=LEGEND_COLUMN_SPACING_EM,
        markerscale=1,
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
        s=SPIKEIN_POINT_SIZE,
        add_equation=True,
        ax=ax,
    )
    ax.set_xlabel(SPIKEIN_X_LABEL)
    ax.set_ylabel(SPIKEIN_Y_LABEL)

    # cnsplots' own legend (with its column-name title) is replaced: the legend
    # above the axes carries the real title.
    add_unframed_legend(ax, columns=frame["Name"].nunique(), title=SPIKEIN_LEGEND_TITLE)


def render_lfc_trajectory_panel(ax: Axes, frame: pd.DataFrame) -> None:
    """Draw one LFC trajectory per gene, coloured by its viability class."""
    for viability in TRAJECTORY_CLASSES:
        shown = frame[frame["FYPOviability"] == viability]
        # One column per gene: a single plot call draws the whole class. The
        # lines carry no label — thousands of them would each become a legend
        # entry — so the legend uses one opaque swatch per class instead.
        curves = shown.pivot_table(index="generation", columns="Systematic ID", values="lfc")
        ax.plot(
            curves.index,
            curves.to_numpy(),
            color=TRAJECTORY_COLORS[viability],
            linewidth=TRAJECTORY_LINE_WIDTH_PT,
            alpha=TRAJECTORY_LINE_ALPHA,
            zorder=TRAJECTORY_ZORDER[viability],
        )

    class_sizes = frame.groupby("FYPOviability")["Systematic ID"].nunique()
    swatches = [
        Line2D(
            [],
            [],
            color=TRAJECTORY_COLORS[viability],
            linewidth=LEGEND_SWATCH_LINE_WIDTH_PT,
            label=f"{TRAJECTORY_LABELS[viability]} (n = {class_sizes[viability]:,})",
        )
        for viability in TRAJECTORY_CLASSES
    ]

    ax.axhline(0.0, color=FURNITURE_COLOR, linestyle="--", linewidth=0.5)
    ax.set_xticks(GENERATIONS, GENERATION_TICK_LABELS)
    ax.set_xlabel(LFC_X_LABEL)
    ax.set_ylabel(LFC_Y_LABEL)
    add_unframed_legend(ax, columns=1, handles=swatches, loc="center", anchor=(0.5, 1.0))


def render_facs_panel(container: Axes, frame: pd.DataFrame) -> None:
    """Draw one mini axes per FACS gene plus the grid's shared labels and legend."""
    genes = sorted(frame["gene"].unique())
    n_rows = facs_gene_rows(frame)
    # One shared, outward-rounded scale lets the ten genes be compared directly.
    y_limits = floor(frame["lfc"].min()), ceil(frame["lfc"].max())
    for index, gene in enumerate(genes):
        row, column = divmod(index, FACS_GRID_COLUMNS)
        ax = container.inset_axes(facs_sub_axes_bounds(row, column, n_rows))
        render_facs_gene_axes(
            ax,
            gene,
            frame.loc[frame["gene"] == gene],
            y_limits,
            show_y_labels=column == 0,
            show_x_labels=row == n_rows - 1,
        )

    add_facs_axis_labels(container, n_rows)
    add_facs_legend(container)


def facs_sub_axes_bounds(row: int, column: int, n_rows: int) -> tuple[float, float, float, float]:
    """Return one mini axes' bounds as fractions of the panel E container."""
    row_from_bottom = n_rows - 1 - row
    x_pt = FACS_LEFT_DECORATION_PT + column * (FACS_MINI_AXES_WIDTH_PT + FACS_COLUMN_GAP_PT)
    y_pt = FACS_BOTTOM_DECORATION_PT + row_from_bottom * (
        FACS_MINI_AXES_HEIGHT_PT + FACS_ROW_GAP_PT
    )
    panel_height_pt = facs_panel_height_pt(n_rows)
    return (
        x_pt / PAGE_CONTENT_WIDTH_PT,
        y_pt / panel_height_pt,
        FACS_MINI_AXES_WIDTH_PT / PAGE_CONTENT_WIDTH_PT,
        FACS_MINI_AXES_HEIGHT_PT / panel_height_pt,
    )


def render_facs_gene_axes(
    ax: Axes,
    gene: str,
    curves: pd.DataFrame,
    y_limits: tuple[int, int],
    show_y_labels: bool,
    show_x_labels: bool,
) -> None:
    """Draw one gene's DIT-HAP and FACS curves under the shared axis conventions."""
    dit_hap = curves[curves["series"] == DIT_HAP_SERIES]
    ax.plot(
        dit_hap["generation"],
        dit_hap["lfc"],
        color=DIT_HAP_COLOR,
        linewidth=FACS_LINE_WIDTH_PT,
        alpha=FACS_LINE_ALPHA,
    )

    for strain in FACS_STRAIN_ORDER:
        facs = curves[curves["series"] == f"{FACS_SERIES_PREFIX}{strain}"]
        if facs.empty:
            continue  # fal1 has only Strain#1; an empty axis would still draw a line
        ax.plot(
            facs["generation"],
            facs["lfc"],
            color=FACS_STRAIN_COLORS[strain],
            linewidth=FACS_LINE_WIDTH_PT,
            alpha=FACS_LINE_ALPHA,
        )

    ax.axhline(0.0, color=FURNITURE_COLOR, linestyle="--", linewidth=0.5)
    ax.set_xlim(*FACS_X_LIMITS)
    ax.set_ylim(*y_limits)
    ax.set_xticks(FACS_X_TICKS)
    lower, upper = y_limits
    ax.set_yticks(range(lower, upper + 1, 2))
    ax.set_title(gene, fontstyle="italic", fontweight="bold")
    # Only the shared labels are turned off; sizes stay at the house 8/7 pt.
    ax.tick_params(
        axis="both",
        labelleft=show_y_labels,
        labelbottom=show_x_labels,
    )


def add_facs_axis_labels(container: Axes, n_rows: int) -> None:
    """Centre one x label and one rotated y label on the whole mini-axes grid."""
    panel_height_pt = facs_panel_height_pt(n_rows)
    grid_height_pt = n_rows * FACS_MINI_AXES_HEIGHT_PT + (n_rows - 1) * FACS_ROW_GAP_PT
    grid_width_pt = (
        FACS_GRID_COLUMNS * FACS_MINI_AXES_WIDTH_PT + (FACS_GRID_COLUMNS - 1) * FACS_COLUMN_GAP_PT
    )
    container.text(
        FACS_Y_LABEL_CENTER_X_PT / PAGE_CONTENT_WIDTH_PT,
        (FACS_BOTTOM_DECORATION_PT + grid_height_pt / 2) / panel_height_pt,
        FACS_Y_LABEL,
        rotation=90,
        ha="center",
        va="center",
    )
    container.text(
        (FACS_LEFT_DECORATION_PT + grid_width_pt / 2) / PAGE_CONTENT_WIDTH_PT,
        FACS_X_LABEL_BOTTOM_PT / panel_height_pt,
        FACS_X_LABEL,
        ha="center",
        va="bottom",
    )


def add_facs_legend(container: Axes) -> None:
    """Draw panel E's three-entry legend as one row above the mini-axes grid."""
    swatches = [
        Line2D(
            [],
            [],
            color=DIT_HAP_COLOR,
            linewidth=FACS_LINE_WIDTH_PT,
            label=DIT_HAP_SERIES,
        ),
        *[
            Line2D(
                [],
                [],
                color=FACS_STRAIN_COLORS[strain],
                linewidth=FACS_LINE_WIDTH_PT,
                label=f"{FACS_SERIES_PREFIX}{strain}",
            )
            for strain in FACS_STRAIN_ORDER
        ],
    ]
    grid_width_pt = (
        FACS_GRID_COLUMNS * FACS_MINI_AXES_WIDTH_PT + (FACS_GRID_COLUMNS - 1) * FACS_COLUMN_GAP_PT
    )
    anchor_x = (FACS_LEFT_DECORATION_PT + grid_width_pt / 2) / PAGE_CONTENT_WIDTH_PT
    add_unframed_legend(
        container,
        columns=len(swatches),
        handles=swatches,
        loc="lower center",
        anchor=(anchor_x, 1.0),
    )


# =============================================================================
# 6. ASSEMBLY (vector overlay and review raster)
# =============================================================================
def page_size_pt(pdf: Path) -> tuple[float, float]:
    """Return the first page's (width, height) in points."""
    page = PdfReader(str(pdf)).pages[0]
    return float(page.mediabox.width), float(page.mediabox.height)


def axes_box_pt(ax: Axes, figure: Figure) -> tuple[float, float, float, float]:
    """Return an axes box as (x, y, width, height) in points from the page's bottom-left."""
    width_in, height_in = figure.get_size_inches()
    width_pt, height_pt = width_in * 72, height_in * 72
    box = ax.get_position()
    return box.x0 * width_pt, box.y0 * height_pt, box.width * width_pt, box.height * height_pt


def assemble_page(
    canvas_pdf: Path,
    schematic_pdf: Path,
    output_pdf: Path,
    box_pt: tuple[float, float, float, float],
) -> None:
    """Overlay the schematic on the canvas, then mount the content atop an A4 page."""
    # The box is in PDF points already: bbox_inches=None keeps the page equal to
    # the figure, so the figure fractions flip straight into the page's origin.
    x_pt, y_pt, box_width_pt, box_height_pt = box_pt
    # The canvas page must belong to a writer before its content can grow: pypdf
    # 7 refuses content replacements on reader-owned pages.
    writer = PdfWriter(clone_from=str(canvas_pdf))
    canvas_page = writer.pages[0]
    schematic_page = PdfReader(str(schematic_pdf)).pages[0]

    source_width = float(schematic_page.mediabox.width)
    source_height = float(schematic_page.mediabox.height)
    scale = min(box_width_pt / source_width, box_height_pt / source_height)

    x_pt += (box_width_pt - source_width * scale) / 2
    y_pt += (box_height_pt - source_height * scale) / 2

    canvas_page.merge_transformed_page(
        schematic_page, Transformation().scale(scale).translate(x_pt, y_pt)
    )

    # A4 page: the canvas keeps its designed size, inset from the page's
    # top-left corner; the rest of the page stays blank.
    content_height_pt = float(canvas_page.mediabox.height)
    bottom_pt = A4_HEIGHT_PT - PAGE_MARGIN_PT - content_height_pt
    page = writer.add_blank_page(width=A4_WIDTH_PT, height=A4_HEIGHT_PT)
    page.merge_transformed_page(canvas_page, Transformation().translate(PAGE_MARGIN_PT, bottom_pt))
    writer.remove_page(0)  # the content now lives on the A4 page
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
def require_columns(
    table: pd.DataFrame, path: Path, required_columns: tuple[str, ...]
) -> pd.DataFrame:
    """Fail once if a table is missing any column the figure reads."""
    missing_columns = set(required_columns) - set(table.columns)
    if missing_columns:
        raise ValueError(f"{path} is missing required columns: {sorted(missing_columns)}")
    return table


def read_tsv(path: Path, required_columns: tuple[str, ...]) -> pd.DataFrame:
    """Read a TSV and validate its columns."""
    return require_columns(pd.read_csv(path, sep="\t"), path, required_columns)


def read_excel(path: Path, required_columns: tuple[str, ...]) -> pd.DataFrame:
    """Read an Excel sheet and validate its columns."""
    return require_columns(pd.read_excel(path), path, required_columns)


def run(config: Figure1Config) -> tuple[float, float]:
    """Build Figure 1 end to end and return the page's (width, height) in points."""
    config.validate()
    apply_house_style()

    coverage = build_coverage_frame(read_tsv(config.lfc_comparison, COVERAGE_COLUMNS))
    spikein = build_spikein_frame(read_tsv(config.spikein, SPIKEIN_COLUMNS))
    gene_lfc = read_tsv(config.gene_lfc, GENE_LFC_COLUMNS)
    trajectory = build_trajectory_frame(gene_lfc)
    facs_comparison = build_facs_comparison_frame(read_excel(config.facs, FACS_COLUMNS), gene_lfc)

    # The schematic page's own box is the panel's box: scale 1.0, actual size.
    schematic_size = page_size_pt(config.schematic)

    mp, schematic_ax = draw_figure(coverage, spikein, trajectory, facs_comparison, schematic_size)
    assert mp.fig is not None
    # cnsplots sets savefig.bbox to "tight", which would crop the page to the
    # drawn artists and shift panel A's box away from the overlay below; keep
    # the designed canvas. The save draws the figure, which is also what settles
    # multipanel's measured layout before the box is read back.
    cns.savefig(config.canvas_pdf, fig=mp.fig, bbox_inches=None)
    assemble_page(
        config.canvas_pdf,
        config.schematic,
        config.output_pdf,
        axes_box_pt(schematic_ax, mp.fig),
    )
    rasterize_review(config.output_pdf, config.output_review)

    logger.success(
        f"Wrote {config.output_pdf} (A4 {A4_WIDTH_PT:.2f} x {A4_HEIGHT_PT:.2f} pt) "
        f"and {config.output_review}"
    )
    return A4_WIDTH_PT, A4_HEIGHT_PT


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Build the composite Figure 1")
    parser.add_argument(
        "--lfc-comparison", type=Path, required=True, help="Combined LFC comparison TSV"
    )
    parser.add_argument("--spikein", type=Path, required=True, help="Spike-in stats TSV")
    parser.add_argument(
        "--gene-lfc", type=Path, required=True, help="Release gene-level LFC TSV (wide)"
    )
    parser.add_argument("--facs", type=Path, required=True, help="FACS depletion-curve workbook")
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
        gene_lfc=args.gene_lfc,
        facs=args.facs,
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
