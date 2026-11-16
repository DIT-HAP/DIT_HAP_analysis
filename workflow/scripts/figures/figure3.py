#!/usr/bin/env python3

"""
Figure 3 — Coherence of annotated gene groups in the DR-DL fitness space
========================================================================

Panel A is the coherence method schematic, drawn outside the pipeline and
supplied as a one-page vector PDF: the observed median pairwise distance in
(DR, DL) space, the permutation null, and the resulting coherence z-score.

Panel B sits to the right of A and asks how significance relates to that
z-score across every annotated term: one point per term in the plane of
coherence z (more negative = tighter than random) against BH-adjusted p on a
-log10 scale. The same two axes carry marginal histograms along the top and the
right, split into all terms (grey) and the coherent cohort (blue), so the shape
of each marginal can be read without leaving the panel. Dashed guides mark the
cohort cut (z < -2, adjusted p <= 0.05); the coherent quadrant is upper left.

Panel C drops the coherent terms onto the fitness axes themselves: one point per
term at the mean DR and mean DL of the genes that were scored for its coherence
test. Marker area is the term's scored size, so a small term and a 200-gene
complex can be told apart. The means are computed here from the release's
per-gene DR/DL table over exactly the members listed in the term's
``scored_member_names`` cell; this is the term-level centre of the same gene set
whose median pairwise distance panel B scores, not the geometric median the
coherence stage publishes (that statistic weights genes differently, and the
panel asks for a plain mean).

All three panels share one row and one 135 pt height: A is a 135 x 135 pt
square, B and C are 150 x 135 pt. A is embedded as vector content, scaled from
its 184.56 pt page box to its 135 pt square, the same pypdf overlay Figure 1
uses for its schematic: matplotlib cannot embed a PDF without rasterizing it.
The review PNG is rasterized from the merged PDF rather than from the canvas.

Both plotted tables come from the de-duplicated representative set (the stage's
paper-facing display layer), so one term contributes one point and the coherent
highlight in B is exactly the set drawn in C.

Every drawing decision lives in this file. The script imports no stage module
from ``workflow/src`` and calls no upstream renderer, so changing Figure 3
cannot change another stage's output.

Input
-----
- results/3a_coherence/<dataset>/dedup/coherence_terms_representatives.tsv
  Columns: source, group_id, group_name, n_scored_members,
  scored_member_names, median_pairwise_distance_z, q_value_coherent.
- <HD_DIT_HAP release>/gene_level/fitting_results.tsv
  Columns: Name, DR, DL (one row per fitted gene).
- resources/figures/figure3a.pdf
  One page whose page box is its content box (184.56 x 184.56 pt); it is
  scaled to the 135 x 135 pt panel square at assembly time.

Output
------
- results/7a_figures/figure3.pdf                content mounted on an A4 page
- results/7a_figures/figure3.review.png         360 dpi review raster
- results/7a_figures/_work/figure3_canvas.pdf   canvas before the A4 mount (debug)

Usage
-----
    python workflow/scripts/figures/figure3.py \\
        --representatives \\
            results/3a_coherence/HD_DIT_HAP/dedup/coherence_terms_representatives.tsv \\
        --fitting-results <HD_DIT_HAP release>/gene_level/fitting_results.tsv \\
        --schematic resources/figures/figure3a.pdf \\
        --output-pdf results/7a_figures/figure3.pdf \\
        --output-review results/7a_figures/figure3.review.png \\
        --canvas results/7a_figures/_work/figure3_canvas.pdf

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-15
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
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from pypdf import PdfReader, PdfWriter, Transformation

# =============================================================================
# 2. CONSTANTS (house style, layout, data vocabulary)
# =============================================================================
# --- House style (self-contained copy; see docs/figures.md) -----------------
HOUSE_PALETTE = "Cell"
FONT_SANS_SERIF = (
    "Arial",
    "Helvetica",
    "Nimbus Sans",
    "Liberation Sans",
    "DejaVu Sans",
)
HOUSE_FIGURE_DPI = 360
FURNITURE_COLOR = cns.GRAY

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

# All three panels share one row and one 135 pt height: A is the square
# schematic box, B and C are 150 pt wide. The wide-but-short B/C keep their
# plotting areas closer to the data's own proportions (a wide volcano, a
# DR-DL cloud with room for DL but no dead vertical space). A's page box is
# scaled to fit its square (vector, via pypdf); B and C lay their insets out in
# their own pt coordinates.
PANEL_A_SIZE_PT = 135.0
PANEL_B_WIDTH_PT = 150.0
PANEL_B_HEIGHT_PT = 135.0
PANEL_B_LEFT_DECORATION_PT = 28.0  # y tick labels + y label
PANEL_B_BOTTOM_DECORATION_PT = 22.0  # x tick labels + x label
PANEL_B_TOP_HISTOGRAM_HEIGHT_PT = 15.0
PANEL_B_RIGHT_HISTOGRAM_WIDTH_PT = 16.0
PANEL_B_HISTOGRAM_GAP_PT = 3.0
PANEL_B_RIGHT_MARGIN_PT = 4.0
PANEL_B_PLOT_WIDTH_PT = (
    PANEL_B_WIDTH_PT
    - PANEL_B_LEFT_DECORATION_PT
    - PANEL_B_RIGHT_HISTOGRAM_WIDTH_PT
    - PANEL_B_HISTOGRAM_GAP_PT
    - PANEL_B_RIGHT_MARGIN_PT
)  # 99

# Panel C is 150 x 135 pt: the plot plus a narrow size-legend band on its right.
PANEL_C_WIDTH_PT = 150.0
PANEL_C_HEIGHT_PT = 135.0
PANEL_C_LEFT_DECORATION_PT = 28.0
PANEL_C_BOTTOM_DECORATION_PT = 22.0
PANEL_C_TOP_MARGIN_PT = 4.0
PANEL_C_PLOT_WIDTH_PT = 92.0
PANEL_C_LEGEND_GAP_PT = 4.0

# Panel B draws the marginal histograms on the same bin grid as the main plot;
# 40 bins keep the count profile smooth at this size without inventing detail.
PANEL_B_BINS = 40

# --- Data vocabulary --------------------------------------------------------
# The de-duplicated representative table is the stage's display layer: one row
# per non-redundant term, with the pooled BH q carried through from combined/.
REPRESENTATIVE_COLUMNS = (
    "source",
    "group_id",
    "group_name",
    "n_scored_members",
    "scored_member_names",
    "median_pairwise_distance_z",
    "q_value_coherent",
)
FITTING_RESULT_COLUMNS = ("Name", "DR", "DL")
SOURCE_COLUMN = "source"
GROUP_ID_COLUMN = "group_id"
GROUP_NAME_COLUMN = "group_name"
SCORED_SIZE_COLUMN = "n_scored_members"
MEMBER_COLUMN = "scored_member_names"
Z_COLUMN = "median_pairwise_distance_z"
Q_COLUMN = "q_value_coherent"
GENE_COLUMN = "Name"
DR_COLUMN = "DR"
DL_COLUMN = "DL"
MEMBER_SEPARATOR = ","

# The cohort definition is config.coherence's; the values are fixed here so the
# figure draws the cut the table was thresholded on.
COHERENT_Q_MAX = 0.05
COHERENT_Z_CUT = -2.0

# --- Panel vocabulary -------------------------------------------------------
COHERENT_COLOR = cns.get_hexcolors_from_apalette([3], HOUSE_PALETTE)[0]  # blue-grey
BACKGROUND_COLOR = "#A3A3A3"  # house grey, the non-coherent cloud
PLOT_TICK_FONT_PT = 6.0
PLOT_LABEL_FONT_PT = 6.5
LEGEND_FONT_PT = 5.0
LEGEND_TITLE_FONT_PT = 5.5
PANEL_B_POINT_SIZE = 4.0
PANEL_B_COHERENT_SIZE = 6.0
PANEL_B_POINT_ALPHA = 0.50
PANEL_B_COHERENT_ALPHA = 0.80
CUTOFF_LINE_WIDTH_PT = 0.5
CUTOFF_LABEL_FONT_PT = 4.5
X_PADDING = 0.03
Y_PADDING = 0.02
SIGNIFICANCE_HEADROOM = 1.03
COHERENCE_TICK_STEP = 2.0
SIGNIFICANCE_TICK_STEP = 1.0
MEAN_DR_TICK_STEP = 0.2
MEAN_DL_TICK_STEP = 1.0
MEAN_DR_PADDING = 0.04
MEAN_DL_HEADROOM = 1.06
PANEL_B_X_LABEL = "Coherence z-score"
PANEL_B_Y_LABEL = r"$-\log_{10}$(BH-adjusted $p$)"
PANEL_C_X_LABEL = "Mean DR"
PANEL_C_Y_LABEL = "Mean DL"
TERM_SIZE_LEGEND_TITLE = "Genes"
TERM_SIZE_LEGEND_VALUES = (5, 20, 50, 100, 200)
TERM_SIZE_SCALE = 0.35  # marker area in pt^2 per scored gene
TERM_POINT_ALPHA = 0.80
TERM_POINT_EDGE_PT = 0.4
TERM_POINT_EDGE_COLOR = "white"


# =============================================================================
# 3. CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class Figure3Config:
    """Resolved input/output paths for Figure 3."""

    representatives: Path
    fitting_results: Path
    schematic: Path
    output_pdf: Path
    output_review: Path
    canvas_pdf: Path

    def validate(self) -> None:
        """Raise if a required input is missing, then create the output dirs."""
        for required in [self.representatives, self.fitting_results, self.schematic]:
            if not required.exists():
                raise ValueError(f"Required input not found: {required}")
        for output in [self.output_pdf, self.output_review, self.canvas_pdf]:
            output.parent.mkdir(parents=True, exist_ok=True)


@dataclass(kw_only=True, slots=True, frozen=True)
class TermCloud:
    """Every representative term's z-score and -log10(q), plus the cohort mask."""

    z: np.ndarray
    significance: np.ndarray
    coherent: np.ndarray


@dataclass(kw_only=True, slots=True, frozen=True)
class CoherentTerm:
    """One coherent term's mean fitness coordinates and scored size."""

    source: str
    group_name: str
    n_scored: int
    mean_dr: float
    mean_dl: float


# =============================================================================
# 4. DATA PREPARATION (pure; schemas are validated at the I/O boundary)
# =============================================================================
def coherent_mask(table: pd.DataFrame) -> np.ndarray:
    """Return the cohort mask defined on the table's q and z columns."""
    is_significant = table[Q_COLUMN].to_numpy(dtype=float) <= COHERENT_Q_MAX
    is_tight = table[Z_COLUMN].to_numpy(dtype=float) < COHERENT_Z_CUT
    return is_significant & is_tight


def build_term_cloud(table: pd.DataFrame) -> TermCloud:
    """Return the panel B arrays for every representative term."""
    z = table[Z_COLUMN].to_numpy(dtype=float)
    q = table[Q_COLUMN].to_numpy(dtype=float)
    if not np.isfinite(z).all() or not np.isfinite(q).all():
        raise ValueError("coherence metrics carry non-finite z or q values")
    if not ((q > 0) & (q <= 1)).all():
        raise ValueError("q_value_coherent must lie in (0, 1] for the -log10 axis")
    return TermCloud(
        z=z,
        significance=-np.log10(q),
        coherent=coherent_mask(table),
    )


def parse_member_names(cell: str) -> tuple[str, ...]:
    """Split a flattened scored-member cell into gene names."""
    return tuple(name.strip() for name in cell.split(MEMBER_SEPARATOR) if name.strip())


def member_frame(table: pd.DataFrame) -> pd.DataFrame:
    """Explode the coherent rows into one (term, member) row per scored gene."""
    coherent = table.loc[coherent_mask(table)].copy()
    if coherent.empty:
        raise ValueError("the representative table has no coherent terms")
    exploded = coherent[
        [SOURCE_COLUMN, GROUP_ID_COLUMN, GROUP_NAME_COLUMN, SCORED_SIZE_COLUMN]
    ].copy()
    exploded[MEMBER_COLUMN] = coherent[MEMBER_COLUMN].map(parse_member_names)
    return exploded.explode(MEMBER_COLUMN).reset_index(drop=True)


def build_coherent_terms(table: pd.DataFrame, fitting: pd.DataFrame) -> tuple[CoherentTerm, ...]:
    """Aggregate each coherent term's scored members to mean DR and mean DL."""
    coherent = table.loc[coherent_mask(table)].copy()
    members = member_frame(table)
    term_keys = [SOURCE_COLUMN, GROUP_ID_COLUMN]

    # The parsed member list must reproduce the size the coherence test used;
    # a mismatch means the flattened cell was truncated or re-formatted.
    expected = coherent.set_index(term_keys, verify_integrity=True)[SCORED_SIZE_COLUMN]
    parsed_counts = members.groupby(term_keys, sort=False).size()
    mismatch = parsed_counts.reindex(expected.index) != expected
    if mismatch.any():
        raise ValueError(
            f"{int(mismatch.sum())} coherent terms list a member count that differs from "
            f"n_scored_members"
        )

    scores = fitting.set_index(GENE_COLUMN)[[DR_COLUMN, DL_COLUMN]]
    if scores.index.duplicated().any():
        raise ValueError("fitting_results.tsv has duplicate gene names")

    joined = members.join(scores, on=MEMBER_COLUMN, how="left")
    missing = joined[DR_COLUMN].isna()
    if missing.any():
        examples = sorted(joined.loc[missing, MEMBER_COLUMN].unique())[:5]
        raise ValueError(f"scored members without a fitted DR/DL: {examples}")

    means = joined.groupby(term_keys, sort=False)[[DR_COLUMN, DL_COLUMN]].mean()
    term_rows = means.reset_index().rename(columns={DR_COLUMN: "mean_dr", DL_COLUMN: "mean_dl"})
    term_rows = term_rows.merge(
        coherent[[SOURCE_COLUMN, GROUP_ID_COLUMN, GROUP_NAME_COLUMN, SCORED_SIZE_COLUMN]],
        on=term_keys,
        how="inner",
        validate="one_to_one",
    )
    return tuple(
        CoherentTerm(
            source=str(row[SOURCE_COLUMN]),
            group_name=str(row[GROUP_NAME_COLUMN]),
            n_scored=int(row[SCORED_SIZE_COLUMN]),
            mean_dr=float(row["mean_dr"]),
            mean_dl=float(row["mean_dl"]),
        )
        for _, row in term_rows.iterrows()
    )


def padded_limits(low: float, high: float, padding: float) -> tuple[float, float]:
    """Pad a data range by a fraction of its span."""
    span = high - low
    return low - padding * span, high + padding * span


def tick_values(low: float, high: float, step: float) -> np.ndarray:
    """Tick the fixed grid that falls inside the plotted range."""
    first = np.ceil(low / step) * step
    last = np.floor(high / step) * step
    return np.arange(first, last + step / 2, step)


def histogram_edges(low: float, high: float, bins: int) -> np.ndarray:
    """Return bin edges over [low, high]."""
    return np.linspace(low, high, bins + 1)


def stacked_counts(
    values: np.ndarray, mask: np.ndarray, edges: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return the background and highlighted counts per bin."""
    all_counts, _ = np.histogram(values, bins=edges)
    coherent_counts, _ = np.histogram(values[mask], bins=edges)
    return all_counts - coherent_counts, coherent_counts


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
    cns.settings.palette_seq = "viridis"
    cns.setup_matplotlib(color_cycle=HOUSE_PALETTE, color_map="viridis")


def draw_figure(cloud: TermCloud, terms: tuple[CoherentTerm, ...]) -> tuple[cns.multipanel, Axes]:
    """Draw panels A, B and C on one row and return the multipanel plus A's axes.

    Returns the multipanel itself, not just its figure: cnsplots keeps the
    measured-layout draw callback as a weak reference, so the caller must hold
    the object until the file is written or the panels stay unmeasured. Panel
    A's axes is returned so the caller can read its settled box for the pypdf
    overlay.
    """
    mp = cns.multipanel(max_width=CANVAS_WIDTH_PT, unit="pt")

    schematic_ax = mp.panel(
        "A",
        width=PANEL_A_SIZE_PT,
        height=PANEL_A_SIZE_PT,
        margin_bottom=0,
    )
    blank_container(schematic_ax)

    container_b = mp.panel(
        "B",
        width=PANEL_B_WIDTH_PT,
        height=PANEL_B_HEIGHT_PT,
        margin_bottom=0,
    )
    blank_container(container_b)
    render_panel_b(container_b, cloud, PANEL_B_HEIGHT_PT)

    container_c = mp.panel(
        "C",
        width=PANEL_C_WIDTH_PT,
        height=PANEL_C_HEIGHT_PT,
        margin_right=0,
        margin_bottom=0,
    )
    blank_container(container_c)
    render_panel_c(container_c, terms, PANEL_C_HEIGHT_PT)
    return mp, schematic_ax


def blank_container(container: Axes) -> None:
    """Clear the ticks, turn the axis off and give the container a pt coordinate system."""
    # Clear the ticks before turning the axis off: multipanel measures y-axis
    # decoration even on an axis-off axes, and placeholder ticks would widen
    # the panel's left reserve.
    container.set_xticks([])
    container.set_yticks([])
    container.set_axis_off()


def render_panel_b(container: Axes, cloud: TermCloud, height_pt: float) -> None:
    """Draw the z-vs-significance cloud and its two marginal histograms."""
    container.set_xlim(0, PANEL_B_WIDTH_PT)
    container.set_ylim(0, height_pt)
    plot_height_pt = (
        height_pt
        - PANEL_B_BOTTOM_DECORATION_PT
        - PANEL_B_TOP_HISTOGRAM_HEIGHT_PT
        - PANEL_B_HISTOGRAM_GAP_PT
    )
    x_limits = padded_limits(float(cloud.z.min()), float(cloud.z.max()), X_PADDING)
    y_limits = padded_limits(0.0, float(cloud.significance.max()), Y_PADDING)
    y_limits = (y_limits[0], y_limits[1] * SIGNIFICANCE_HEADROOM)

    plot_ax = container.inset_axes(
        (
            PANEL_B_LEFT_DECORATION_PT / PANEL_B_WIDTH_PT,
            PANEL_B_BOTTOM_DECORATION_PT / height_pt,
            PANEL_B_PLOT_WIDTH_PT / PANEL_B_WIDTH_PT,
            plot_height_pt / height_pt,
        )
    )
    top_ax = container.inset_axes(
        (
            PANEL_B_LEFT_DECORATION_PT / PANEL_B_WIDTH_PT,
            (PANEL_B_BOTTOM_DECORATION_PT + plot_height_pt + PANEL_B_HISTOGRAM_GAP_PT) / height_pt,
            PANEL_B_PLOT_WIDTH_PT / PANEL_B_WIDTH_PT,
            PANEL_B_TOP_HISTOGRAM_HEIGHT_PT / height_pt,
        )
    )
    right_ax = container.inset_axes(
        (
            (PANEL_B_LEFT_DECORATION_PT + PANEL_B_PLOT_WIDTH_PT + PANEL_B_HISTOGRAM_GAP_PT)
            / PANEL_B_WIDTH_PT,
            PANEL_B_BOTTOM_DECORATION_PT / height_pt,
            PANEL_B_RIGHT_HISTOGRAM_WIDTH_PT / PANEL_B_WIDTH_PT,
            plot_height_pt / height_pt,
        )
    )
    render_panel_b_plot(plot_ax, cloud, x_limits, y_limits)
    render_panel_b_marginals(top_ax, right_ax, plot_ax, cloud, x_limits, y_limits)


def render_panel_b_plot(
    ax: Axes,
    cloud: TermCloud,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
) -> None:
    """Scatter every term and mark the cohort cut on both axes."""
    is_background = ~cloud.coherent
    ax.scatter(
        cloud.z[is_background],
        cloud.significance[is_background],
        s=PANEL_B_POINT_SIZE,
        color=BACKGROUND_COLOR,
        alpha=PANEL_B_POINT_ALPHA,
        linewidths=0,
    )
    ax.scatter(
        cloud.z[cloud.coherent],
        cloud.significance[cloud.coherent],
        s=PANEL_B_COHERENT_SIZE,
        color=COHERENT_COLOR,
        alpha=PANEL_B_COHERENT_ALPHA,
        linewidths=0,
    )
    ax.axvline(
        COHERENT_Z_CUT,
        color=FURNITURE_COLOR,
        linestyle="--",
        linewidth=CUTOFF_LINE_WIDTH_PT,
        zorder=0,
    )
    ax.axhline(
        -np.log10(COHERENT_Q_MAX),
        color=FURNITURE_COLOR,
        linestyle="--",
        linewidth=CUTOFF_LINE_WIDTH_PT,
        zorder=0,
    )
    ax.text(
        x_limits[1],
        -np.log10(COHERENT_Q_MAX),
        f"$q = {COHERENT_Q_MAX:g}$",
        ha="right",
        va="bottom",
        fontsize=CUTOFF_LABEL_FONT_PT,
        color=FURNITURE_COLOR,
    )
    ax.text(
        COHERENT_Z_CUT,
        y_limits[1],
        f"$z = {COHERENT_Z_CUT:g}$",
        ha="left",
        va="top",
        rotation=90,
        fontsize=CUTOFF_LABEL_FONT_PT,
        color=FURNITURE_COLOR,
    )
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_xticks(tick_values(*x_limits, COHERENCE_TICK_STEP))
    ax.set_yticks(tick_values(y_limits[0], y_limits[1], SIGNIFICANCE_TICK_STEP))
    ax.set_xlabel(PANEL_B_X_LABEL, fontsize=PLOT_LABEL_FONT_PT)
    ax.set_ylabel(PANEL_B_Y_LABEL, fontsize=PLOT_LABEL_FONT_PT)
    ax.tick_params(axis="both", labelsize=PLOT_TICK_FONT_PT)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    add_panel_b_legend(ax, cloud)


def add_panel_b_legend(ax: Axes, cloud: TermCloud) -> None:
    """Name the two point styles in the empty upper-right corner."""
    handles = [
        Line2D([], [], marker="o", linestyle="", markersize=3.0, color=BACKGROUND_COLOR),
        Line2D([], [], marker="o", linestyle="", markersize=3.5, color=COHERENT_COLOR),
    ]
    labels = [
        f"All terms (n = {len(cloud.z):,})",
        f"Coherent (n = {int(cloud.coherent.sum()):,})",
    ]
    ax.legend(
        handles,
        labels,
        loc="upper right",
        frameon=False,
        fontsize=LEGEND_FONT_PT,
        handlelength=0.8,
        handletextpad=0.4,
        labelspacing=0.3,
        borderaxespad=0.2,
    )


def render_panel_b_marginals(
    top_ax: Axes,
    right_ax: Axes,
    plot_ax: Axes,
    cloud: TermCloud,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
) -> None:
    """Draw the z and significance distributions along their own axes."""
    z_edges = histogram_edges(x_limits[0], x_limits[1], PANEL_B_BINS)
    significance_edges = histogram_edges(y_limits[0], y_limits[1], PANEL_B_BINS)
    z_background, z_coherent = stacked_counts(cloud.z, cloud.coherent, z_edges)
    q_background, q_coherent = stacked_counts(
        cloud.significance, cloud.coherent, significance_edges
    )

    # Empty bins are skipped: a zero-height bar still renders as a hairline at
    # the baseline, which would draw a false rule across the histogram.
    z_occupied = (z_background + z_coherent) > 0
    z_highlighted = z_coherent > 0
    top_ax.bar(
        z_edges[:-1][z_occupied],
        z_background[z_occupied],
        width=np.diff(z_edges)[z_occupied],
        align="edge",
        color=BACKGROUND_COLOR,
        linewidth=0,
    )
    top_ax.bar(
        z_edges[:-1][z_highlighted],
        z_coherent[z_highlighted],
        bottom=z_background[z_highlighted],
        width=np.diff(z_edges)[z_highlighted],
        align="edge",
        color=COHERENT_COLOR,
        linewidth=0,
    )
    q_occupied = (q_background + q_coherent) > 0
    q_highlighted = q_coherent > 0
    right_ax.barh(
        significance_edges[:-1][q_occupied],
        q_background[q_occupied],
        height=np.diff(significance_edges)[q_occupied],
        align="edge",
        color=BACKGROUND_COLOR,
        linewidth=0,
    )
    right_ax.barh(
        significance_edges[:-1][q_highlighted],
        q_coherent[q_highlighted],
        left=q_background[q_highlighted],
        height=np.diff(significance_edges)[q_highlighted],
        align="edge",
        color=COHERENT_COLOR,
        linewidth=0,
    )

    top_ax.axvline(
        COHERENT_Z_CUT,
        color=FURNITURE_COLOR,
        linestyle="--",
        linewidth=CUTOFF_LINE_WIDTH_PT,
    )
    right_ax.axhline(
        -np.log10(COHERENT_Q_MAX),
        color=FURNITURE_COLOR,
        linestyle="--",
        linewidth=CUTOFF_LINE_WIDTH_PT,
    )
    top_ax.sharex(plot_ax)
    right_ax.sharey(plot_ax)
    top_ax.set_ylim(0, float((z_background + z_coherent).max() * 1.08))
    right_ax.set_xlim(0, float((q_background + q_coherent).max() * 1.08))
    hide_frame(top_ax)
    hide_frame(right_ax)


def hide_frame(ax: Axes) -> None:
    """Remove every frame element from a marginal axes, sharing axes included."""
    ax.set_axis_off()
    ax.set_frame_on(False)
    ax.patch.set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def render_panel_c(container: Axes, terms: tuple[CoherentTerm, ...], height_pt: float) -> None:
    """Draw one mean-DR/mean-DL bubble per coherent term."""
    container.set_xlim(0, PANEL_C_WIDTH_PT)
    container.set_ylim(0, height_pt)
    plot_height_pt = height_pt - PANEL_C_BOTTOM_DECORATION_PT - PANEL_C_TOP_MARGIN_PT
    ax = container.inset_axes(
        (
            PANEL_C_LEFT_DECORATION_PT / PANEL_C_WIDTH_PT,
            PANEL_C_BOTTOM_DECORATION_PT / height_pt,
            PANEL_C_PLOT_WIDTH_PT / PANEL_C_WIDTH_PT,
            plot_height_pt / height_pt,
        )
    )
    mean_dr = np.array([term.mean_dr for term in terms])
    mean_dl = np.array([term.mean_dl for term in terms])
    sizes = np.array([term.n_scored for term in terms], dtype=float)
    x_limits = padded_limits(float(mean_dr.min()), float(mean_dr.max()), MEAN_DR_PADDING)
    y_limits = (0.0, float(mean_dl.max()) * MEAN_DL_HEADROOM)

    ax.scatter(
        mean_dr,
        mean_dl,
        s=sizes * TERM_SIZE_SCALE,
        color=COHERENT_COLOR,
        alpha=TERM_POINT_ALPHA,
        edgecolors=TERM_POINT_EDGE_COLOR,
        linewidths=TERM_POINT_EDGE_PT,
    )
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_xticks(tick_values(*x_limits, MEAN_DR_TICK_STEP))
    ax.set_yticks(tick_values(y_limits[0], y_limits[1], MEAN_DL_TICK_STEP))
    ax.set_xlabel(PANEL_C_X_LABEL, fontsize=PLOT_LABEL_FONT_PT)
    ax.set_ylabel(PANEL_C_Y_LABEL, fontsize=PLOT_LABEL_FONT_PT)
    ax.tick_params(axis="both", labelsize=PLOT_TICK_FONT_PT)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    add_panel_c_size_legend(ax)


def add_panel_c_size_legend(ax: Axes) -> None:
    """Draw the marker-area key to the right of the panel C plot."""
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            markerfacecolor=COHERENT_COLOR,
            markeredgecolor=TERM_POINT_EDGE_COLOR,
            markeredgewidth=TERM_POINT_EDGE_PT,
            markersize=np.sqrt(size * TERM_SIZE_SCALE),
        )
        for size in TERM_SIZE_LEGEND_VALUES
    ]
    labels = [f"{size:,}" for size in TERM_SIZE_LEGEND_VALUES]
    legend_x = (PANEL_C_PLOT_WIDTH_PT + PANEL_C_LEGEND_GAP_PT) / PANEL_C_PLOT_WIDTH_PT
    ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(legend_x, 1.0),
        frameon=False,
        fontsize=LEGEND_FONT_PT,
        title=TERM_SIZE_LEGEND_TITLE,
        title_fontsize=LEGEND_TITLE_FONT_PT,
        handletextpad=0.5,
        labelspacing=0.6,
        borderaxespad=0.0,
    )


# =============================================================================
# 6. ASSEMBLY
# =============================================================================
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
def read_tsv_columns(path: Path, required: tuple[str, ...]) -> pd.DataFrame:
    """Read only the required TSV columns and fail once if any is missing."""
    table = pd.read_csv(path, sep="\t", usecols=lambda column: column in required)
    missing_columns = set(required) - set(table.columns)
    if missing_columns:
        raise ValueError(f"{path} is missing required columns: {sorted(missing_columns)}")
    return table


def run(config: Figure3Config) -> tuple[float, float]:
    """Build Figure 3 end to end and return the page's (width, height) in points."""
    config.validate()
    apply_house_style()

    representatives = read_tsv_columns(config.representatives, REPRESENTATIVE_COLUMNS)
    fitting = read_tsv_columns(config.fitting_results, FITTING_RESULT_COLUMNS)
    cloud = build_term_cloud(representatives)
    terms = build_coherent_terms(representatives, fitting)
    logger.info(
        f"Panel B: {len(cloud.z):,} terms ({int(cloud.coherent.sum()):,} coherent); "
        f"panel C: {len(terms):,} coherent terms"
    )

    mp, schematic_ax = draw_figure(cloud, terms)
    assert mp.fig is not None
    # cnsplots sets savefig.bbox to "tight", which would crop the page to the
    # drawn artists and shift panel A's box away from the overlay below; keep
    # the designed canvas. The save draws the figure, which is also what
    # settles multipanel's measured layout before the box is read back.
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
    parser = argparse.ArgumentParser(description="Build the composite Figure 3")
    parser.add_argument(
        "--representatives",
        type=Path,
        required=True,
        help="De-duplicated representative terms TSV",
    )
    parser.add_argument(
        "--fitting-results",
        type=Path,
        required=True,
        help="Release per-gene fitting results (Name + DR + DL)",
    )
    parser.add_argument("--schematic", type=Path, required=True, help="Panel A schematic PDF")
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
    config = Figure3Config(
        representatives=args.representatives,
        fitting_results=args.fitting_results,
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
