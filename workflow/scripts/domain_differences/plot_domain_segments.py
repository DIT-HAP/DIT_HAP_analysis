#!/usr/bin/env python3

"""
Plot Domain Segments
====================

Stage 4a, step 5: the stage's one figure — a gallery with one gene per row, six
genes per page, one page per group of the hand annotation.

A row is read like this:

- blue dot = a Forward insertion, red dot = a Reverse insertion, thin grey tick
  joining the pair; the open diamond on top is their mean, and that mean is the
  value the segmentation actually ran on.
- the coloured horizontal bars are the segments, each drawn at its own ``mean_DR``;
  dashed verticals mark the boundaries.
- the grey band is the merge threshold expressed in DR units, ``lambda_dr`` — the
  paper's detail threshold ``C*sqrt(2(1+delta)lnT)*sigma_gene + sigma_add``
  divided by ``c(nL, nR)``, taken at the median surviving boundary — centred on
  the gene's own level: a step that fits inside it is noise to this rule, so a
  segment only means something if its bar sits outside the band.
- green shading marks the hand-curated non-essential region, where there is one.

Above the DR profile sit three AFdb domain tracks on the same residue axis:
PAE (AlphaFold PAE domains, disordered calls omitted), TED (pLDDT segments;
``medium`` drawn lighter than ``high``) and ECOD (structural domains). Each
domain gets its own colour — the house palette, then seeded random hues past its
end — and a domain split into several parts is drawn as its parts, so a gap
stays a gap. A track a gene has no call for, or that the upstream table does not
cover at all, reads ``n/a``. Colouring is keyed to the domain, not to a residue
range: the same ``segment_id`` keeps one colour across all its parts.

Every row carries that key itself, in one column beside the panel, and a two-line
caption under it: ``name segs=3 n_sites=8 sigma=0.24 lambda=0.47`` over the first
few segment ranges and ``max|Δ|``. Every panel also titles itself with the gene,
its ``FYPOviability`` call and its gRNA-level DR, and labels both its axes — a
reader who has scrolled to one gene needs nothing from anywhere else on the page.
Only the group heading is shared, and it sits at the top of the page.

Genes are grouped C -> N -> N+C -> unannotated and sorted by
``max_boundary_contrast`` inside each group, so the hand annotation — the only
external check the stage has — is what a reader sees first.

Input
-----
- domain_summary.tsv, domain_segments.tsv — the stage's tables
- domain_hand_regions.tsv — the hand annotation
- _work/insertions.parquet — the raw Forward/Reverse points
- gene_annotation_reference.protein.parquet (1c) — the panel headings' annotation
- resources/external/afdb_domains/{version}/afdb_domain_segments.tsv — the
  upstream AFdb PAE/TED/ECOD segments, positioned on the same PomBase peptide

Output
------
- domain_segments.pdf — the gallery
- domain_segments.review.png — page 1, for a quick look
- domain_segments/<group>_<n>_genes.pdf — the same pages split by gene group, the
  group's gene count in the file name; each also gets its page-1 review PNG. Every
  page heading gives the group, the gene range on that page and the group total
  (e.g. ``C terminus — 1-6 / 102``).

Usage
-----
    python workflow/scripts/domain_differences/plot_domain_segments.py \\
        --segments results/4a_domain_differences/HD_DIT_HAP/domain_segments.tsv \\
        --summary results/4a_domain_differences/HD_DIT_HAP/domain_summary.tsv \\
        --hand-regions results/4a_domain_differences/HD_DIT_HAP/domain_hand_regions.tsv \\
        --insertions results/4a_domain_differences/HD_DIT_HAP/_work/insertions.parquet \\
        --annotation
        results/1c_annotation/2026-06-01/2026-08-11/gene_annotation_reference.protein.parquet \\
        --afdb-segments
        resources/external/afdb_domains/2026-11-16/afdb_domain_segments.tsv \\
        --output results/4a_domain_differences/HD_DIT_HAP/domain_segments.pdf

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-10-14
Version:  1.1.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
from collections.abc import Iterator
import colorsys
from dataclasses import dataclass
from pathlib import Path
import random
import sys

# 2. Third-party Imports
from loguru import logger

# 3. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
import cnsplots as cns
from matplotlib.axes import Axes
from matplotlib.colors import to_hex
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from domain_differences.core import (
    TERM_ORDER,
    UNANNOTATED,
    collapse_sites,
    hand_region_map,
    read_afdb_domains,
    read_gene_context,
)
from figures import (
    FIT_CAPTION_FONT_OFFSET_PT,
    FURNITURE_COLOR,
    HOUSE_SCATTER_S,
    PANEL_SIZES_PX,
    PanelShape,
    apply_house_style,
    cluster_color,
    fit_panels,
    grid_axes,
    house_palette_length,
    save_pdf_pages,
)
from io_table import read_parquet
from logging_setup import setup_logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
ROWS_PER_PAGE = 6
GENE_GROUP_LABELS = {**{T: f"{T} terminus" for T in TERM_ORDER}, UNANNOTATED: "no hand annotation"}
# File-name slugs for the per-group PDFs, which also carry their gene count.
GENE_GROUP_SLUGS = {
    "C": "c_terminus",
    "N": "n_terminus",
    "N+C": "n_and_c_termini",
    UNANNOTATED: "unannotated",
}
# Residue ranges spelled out under a panel before they are abbreviated. Three are
# what fits on one caption line without pushing max|Δ| off the panel; the count
# itself is on the line above, so a truncated list still says how many there are.
SEGMENTS_SHOWN = 3
# Width of the legend column every panel carries on its right, reserved on the
# page BEFORE fit_panels measures. fit_panels corrects a grid whose cells are the
# wrong size by growing the page; a panel has a fixed box aspect, so a first-pass
# width error arrives as an equally large height error, and on a 6-row page that
# means growing the height 6 rows deep — space fit_panels never takes back, left
# blank between the rows. Reserving the column keeps that first measurement honest.
LEGEND_COLUMN_IN = 1.45
# Page fraction reserved above the grid for the group heading. Every panel titles
# itself now, so the group name has to live at the page level, and tight_layout
# only reserves room for a suptitle on the pass that added it — fit_panels then
# grows the page and lays out again, which on its own left the group heading
# touching the first panel's title.
SUPTITLE_STRIP = 0.02

# DR spans [-2.09, 0.17] over the whole release; one fixed range for every panel is
# what lets the threshold bands be compared between rows, and it covers the data
# instead of clipping its tail.
Y_RANGE = (-2.2, 0.4)
# The two caption lines take about a sixth of a 100 px panel. Reserving that
# fraction at the foot — by widening the y range downward, not by shrinking the
# panel — is what keeps the text off the points: the release's lowest insertion
# sits at -2.09, INSIDE the window, so a caption drawn over the window covers it.
CAPTION_BAND = 0.17

# The panel is a WIDE box plus a strip on top. The 30 pt strip holds the AFdb
# domain calls in three lanes; the DR profile keeps the 100 pt it had, so the
# strip costs page height (130/100) rather than data resolution.
PANEL_SHAPE = PanelShape.WIDE_TALL
TRACK_BAND_PT = 30.0
TRACK_LANES = ("PAE", "TED", "ECOD")
TRACK_BAND = TRACK_BAND_PT / PANEL_SIZES_PX[PANEL_SHAPE][1]
TRACK_LANE_HEIGHT = TRACK_BAND / len(TRACK_LANES)
TRACK_BAR_LW = 0.7 * TRACK_LANE_HEIGHT * PANEL_SIZES_PX[PANEL_SHAPE][1]
# Probed once: house_palette_length() raises its way up to the length, and the
# tracks would otherwise call it per gene and lane.
PALETTE_LENGTH = house_palette_length()
# TED calls its segments high/medium confidence; the medium ones keep the domain's
# hue and are drawn lighter rather than in a different colour.
TED_CLASS_ALPHA = {"high": 1.0, "medium": 0.45}

# The y window: the data range plus its caption band keeps its old DR-unit span,
# and the strip above it is added as the remaining axes fraction.
_DATA_SPAN = (Y_RANGE[1] - Y_RANGE[0]) / (1 - CAPTION_BAND)
_Y_BOTTOM = Y_RANGE[0] - CAPTION_BAND / (1 - CAPTION_BAND) * (Y_RANGE[1] - Y_RANGE[0])
Y_TOP = _Y_BOTTOM + _DATA_SPAN / (1 - TRACK_BAND)
# Ticks stay inside the data window; without an explicit list matplotlib would
# place them across the track strip too.
Y_TICKS = np.arange(np.floor(Y_RANGE[0] * 2) / 2, Y_RANGE[1], 0.5)

# One colour per meaning across the whole figure. The Forward/Reverse pair and the
# hand annotation come from the house palette; the site mean is deliberately
# neutral so the two experimental channels stay the loudest thing in a panel.
FWD_COLOR = cns.BLUE
REV_COLOR = cns.RED
SEGMENT_COLOR = cns.PURPLE
HAND_COLOR = cns.GREEN
SITE_MEAN_COLOR = "0.25"
LAMBDA_COLOR = FURNITURE_COLOR


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Runtime parameters for the gallery."""

    segments: Path
    summary: Path
    hand_regions: Path
    insertions: Path
    annotation: Path
    afdb_segments: Path
    output: Path


@dataclass(kw_only=True, slots=True, frozen=True)
class Gallery:
    """Everything one render needs, grouped and indexed for lookup by gene."""

    genes: pd.DataFrame
    segments: dict
    sites: dict
    raw: dict
    hand: dict
    domains: dict
    viability: dict
    grna_dr: dict


# =============================================================================
# CORE LOGIC — data for the panels
# =============================================================================
def load_gallery(config: PlotConfig) -> Gallery:
    """Read the stage's tables and index them by gene, in gallery order."""
    segments = pd.read_csv(config.segments, sep="\t")
    hand_table = pd.read_csv(config.hand_regions, sep="\t")
    hand = hand_region_map(hand_table)

    insertions = read_parquet(config.insertions)
    # list() before dict(): pandas 2.3 exposes GroupBy.keys as the grouping
    # column name, so dict() mistakes the groupby for a mapping and calls that
    # string.
    raw = dict(list(insertions.groupby("Systematic ID", sort=False)))
    sites = {sid: collapse_sites(frame) for sid, frame in raw.items()}
    viability, grna_dr = read_gene_context(config.annotation)

    # The term a gene is grouped under. No gene is annotated in two sheets, but if
    # that ever changed the first sheet in TERM_ORDER wins.
    term = (
        hand_table.dropna(subset=["sid"]).drop_duplicates("sid").set_index("sid")["term"].to_dict()
    )
    rank = {t: i for i, t in enumerate((*TERM_ORDER, UNANNOTATED))}
    summary = pd.read_csv(config.summary, sep="\t")
    summary["term"] = [term.get(sid, UNANNOTATED) for sid in summary["Systematic ID"]]
    summary["_rank"] = summary["term"].map(rank).astype(int)
    genes = (
        summary.sort_values(["_rank", "max_boundary_contrast"], ascending=[True, False])
        .drop(columns="_rank")
        .reset_index(drop=True)
    )
    protein_lengths = {
        sid: int(n) for sid, n in zip(summary["Systematic ID"], summary["protein_len"], strict=True)
    }
    domains = read_afdb_domains(config.afdb_segments, protein_lengths)

    return Gallery(
        genes=genes,
        segments=dict(list(segments.groupby("Systematic ID", sort=False))),
        sites=sites,
        raw=raw,
        hand=hand,
        domains=domains,
        viability=viability,
        grna_dr=grna_dr,
    )


# =============================================================================
# CORE LOGIC — one panel
# =============================================================================
def draw_insertions(ax: Axes, raw: pd.DataFrame, sites: pd.DataFrame) -> None:
    """Raw Forward/Reverse points of one gene plus their per-site mean.

    The typical Forward-vs-Reverse gap is only 0.045 DR, a couple of pixels at this
    page scale, so a pair drawn as two dots looks like one point. The thin vertical
    tick joining them makes every site read as a pair, and its length IS the
    disagreement. The mean is drawn as an OPEN diamond so it never hides the pair.
    """
    forward = raw["Insertion_direction"].eq("Forward").to_numpy()
    residues = raw["Residue_affected"].astype(int).to_numpy()
    ax.vlines(
        sites["Residue_affected"].astype(int),
        sites["DR_rev"].fillna(sites["DR"]),
        sites["DR_fwd"].fillna(sites["DR"]),
        color="0.6",
        zorder=3,
    )
    ax.scatter(
        residues[forward],
        raw["DR"].to_numpy()[forward],
        s=HOUSE_SCATTER_S,
        color=FWD_COLOR,
        zorder=4,
    )
    ax.scatter(
        residues[~forward],
        raw["DR"].to_numpy()[~forward],
        s=HOUSE_SCATTER_S,
        color=REV_COLOR,
        zorder=4,
    )
    ax.scatter(
        sites["Residue_affected"].astype(int),
        sites["DR"],
        s=HOUSE_SCATTER_S,
        marker="D",
        facecolors="none",
        edgecolors=SITE_MEAN_COLOR,
        zorder=5,
    )


def domain_colors(intervals: pd.DataFrame, rng: random.Random) -> dict[str, str]:
    """One colour per structure domain: the house palette, then seeded random hues."""
    order = intervals.groupby("segment_id")["start"].min().sort_values().index
    in_palette = min(len(order), PALETTE_LENGTH)
    colors = [cluster_color(k) for k in range(1, in_palette + 1)]
    colors += [
        to_hex(colorsys.hsv_to_rgb(rng.random(), 0.55, 0.75))
        for _ in range(len(order) - in_palette)
    ]
    return dict(zip(order, colors, strict=True))


def draw_domain_tracks(
    ax: Axes, domains: pd.DataFrame | None, caption_size: float, sid: str
) -> None:
    """PAE/TED/ECOD calls in three lanes above the DR profile, on the panel's x axis.

    A domain split into several parts is drawn as its parts, so a gap stays a gap.
    PAE's disordered calls are left out: the lane is a domain map, and a gap then
    reads as low structural confidence rather than as a missing call.
    """
    # Seeded per gene so a rebuild paints the same random colours past the palette.
    rng = random.Random(sid)
    for index, source in enumerate(TRACK_LANES):
        # index 0 is the top lane: axes fraction counts up from the bottom
        lane = 1 - TRACK_LANE_HEIGHT * (index + 0.5)
        ax.text(
            -0.02,
            lane,
            source,
            transform=ax.transAxes,
            ha="right",
            va="center",
            color=SITE_MEAN_COLOR,
            fontsize=caption_size,
        )
        intervals = None if domains is None else domains.loc[domains["source"].eq(source)]
        if intervals is not None and source == "PAE":
            intervals = intervals.loc[intervals["segment_class"].eq("domain")]
        if intervals is None or intervals.empty:
            ax.text(
                0.01,
                lane,
                "n/a",
                transform=ax.transAxes,
                ha="left",
                va="center",
                color=FURNITURE_COLOR,
                fontsize=caption_size,
            )
            continue

        colors = domain_colors(intervals, rng)
        for _, interval in intervals.iterrows():
            start = float(interval["start"])
            end = float(interval["end"])
            color = colors[interval["segment_id"]]
            alpha = TED_CLASS_ALPHA.get(interval["segment_class"], 1.0) if source == "TED" else 1.0
            if start == end:
                # a one-residue call has zero width in data units; a tick is its footprint
                ax.plot(
                    [start],
                    [lane],
                    marker="|",
                    ls="",
                    markersize=TRACK_BAR_LW * 0.6,
                    transform=ax.get_xaxis_transform(),
                    color=color,
                    alpha=alpha,
                    zorder=4,
                )
            else:
                ax.plot(
                    [start, end],
                    [lane, lane],
                    lw=TRACK_BAR_LW,
                    solid_capstyle="butt",
                    transform=ax.get_xaxis_transform(),
                    color=color,
                    alpha=alpha,
                    zorder=4,
                )


def panel_title(name: str, viability: str | None, grna_dr: float | None) -> str:
    """The panel's heading: what the gene is, two lines because one will not fit.

    A name, a viability call and a gRNA DR run to 55 characters at the long end
    (``SPBC1861.05  —  depends_on_conditions  —  gRNA DR=-1.04``) — more than a
    150 pt panel holds at any legible size — so the name and the call take one
    line and the DR the next.

    Both annotations come from the reference the rest of the project reads; a
    gene missing from it shows ``—`` rather than being dropped, since the gallery
    is built from the segmentation, not from the annotation.
    """
    call = viability or "—"
    # `grna_dr == grna_dr` is False only for NaN, which read_gene_context already
    # filters out — kept because a NaN that printed as "gRNA DR=+nan" would look
    # like a value rather than a missing one
    dr = f"gRNA DR={grna_dr:+.2f}" if grna_dr is not None and grna_dr == grna_dr else "gRNA DR=—"
    return f"{name}  —  {call}\n{dr}"


def draw_panel(
    ax: Axes,
    row: pd.Series,
    segments: pd.DataFrame,
    raw: pd.DataFrame,
    sites: pd.DataFrame,
    hand: list[tuple[int, int]],
    domains: pd.DataFrame | None = None,
    viability: str | None = None,
    grna_dr: float | None = None,
) -> None:
    """One gene: heading, threshold band, hand shading, the two channels, segment levels,
    domain tracks, caption."""
    # Registered fit-critical size: the gallery's widest caption runs 149 pt at
    # 6.5 pt into a 150 pt panel and crosses the spine; 6.0 leaves 9 pt of margin.
    caption_size = plt.rcParams["font.size"] + FIT_CAPTION_FONT_OFFSET_PT
    ax.set_title(panel_title(str(row["name"]), viability, grna_dr))
    ax.axhline(0, color=LAMBDA_COLOR, zorder=0)

    lam = float(row["lambda_dr"])
    level = float(sites["DR"].mean())
    if lam > 0:
        ax.axhspan(level - lam / 2, level + lam / 2, color=LAMBDA_COLOR, alpha=0.13, zorder=1)
        for edge in (level - lam / 2, level + lam / 2):
            ax.axhline(edge, color=LAMBDA_COLOR, ls=":", zorder=2)
        ax.text(
            0.98,
            level + lam / 2,
            f"λ_DR={lam:.2f}",
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="bottom",
            color=SITE_MEAN_COLOR,
            fontsize=caption_size,
        )

    for start, end in hand:
        ax.axvspan(start, end, color=HAND_COLOR, alpha=0.18, zorder=1, ymax=1 - TRACK_BAND)

    draw_insertions(ax, raw, sites)
    draw_domain_tracks(ax, domains, caption_size, str(row["Systematic ID"]))
    for _, segment in segments.iterrows():
        ax.hlines(
            float(segment["mean_DR"]),
            int(segment["start_residue"]),
            int(segment["end_residue"]),
            color=SEGMENT_COLOR,
            zorder=6,
        )
    for _, segment in segments.iloc[1:].iterrows():
        ax.axvline(int(segment["start_residue"]), color="0.55", ls="--", zorder=5)

    # The x range is the whole protein, not just the measured part: the tracks
    # annotate it end to end, so a window cut to the insertions would clip them.
    # The unmeasured flanks then read as the absence of points, which is what they are.
    ax.set_xlim(0.5, float(row["protein_len"]) + 0.5)
    # see CAPTION_BAND / TRACK_BAND: the data window keeps its old footprint at the
    # foot, the track strip owns the fraction above it
    ax.set_ylim(_Y_BOTTOM, Y_TOP)
    ax.set_yticks(Y_TICKS)

    ranges = "|".join(
        f"{int(s['start_residue'])}-{int(s['end_residue'])}"
        for _, s in segments.iloc[:SEGMENTS_SHOWN].iterrows()
    )
    if len(segments) > SEGMENTS_SHOWN:
        ranges += f"|+{len(segments) - SEGMENTS_SHOWN}"
    ax.text(
        0.02,
        0.04,
        f"{row['name']}  segs={int(row['n_segments'])}  n_sites={int(row['n_sites'])}"
        f"  \u03c3={float(row['sigma_gene']):.2f}  λ_DR={lam:.2f}\n"
        f"{ranges}  max|Δ|={float(row['max_boundary_contrast']):.2f}",
        transform=ax.transAxes,
        va="bottom",
        linespacing=1.3,
        fontsize=caption_size,
    )

    # every row carries the key, in one column beside it: a reader who scrolled to
    # a single gene should not have to find the page's legend to know which dots
    # are Forward. Outside the axes, so it never sits on the data.
    ax.legend(
        handles=legend_handles(),
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        ncol=1,
        frameon=False,
        fontsize=caption_size,
        handletextpad=0.6,
        borderaxespad=0.0,
        labelspacing=0.8,
    )


def legend_handles() -> list[Line2D]:
    """One key per visual channel, with the site mean open like it is in the panels."""
    return [
        Line2D([], [], marker="o", ls="", color=FWD_COLOR, label="insertion Forward"),
        Line2D([], [], marker="o", ls="", color=REV_COLOR, label="insertion Reverse"),
        Line2D(
            [],
            [],
            marker="D",
            ls="",
            markerfacecolor="none",
            markeredgecolor=SITE_MEAN_COLOR,
            label="site mean (segmented on)",
        ),
        Line2D([], [], color=SEGMENT_COLOR, lw=1.6, label="segment mean DR"),
        Line2D([], [], color=LAMBDA_COLOR, lw=6, alpha=0.35, label="λ_DR at surviving boundaries"),
        Line2D([], [], color=HAND_COLOR, lw=6, alpha=0.35, label="hand-curated non-essential"),
    ]


# =============================================================================
# CORE LOGIC — the gallery
# =============================================================================
def draw_page(page: pd.DataFrame, gallery: Gallery, start: int, end: int, total: int) -> None:
    """Draw one page of up to ROWS_PER_PAGE genes into a freshly created figure."""
    # empty labels: a row is a named gene, not a lettered panel, and its own title
    # already carries the identity
    axes = grid_axes(ROWS_PER_PAGE, 1, labels=[], shape=PANEL_SHAPE)
    for ax, (_, row) in zip(axes, page.iterrows(), strict=False):
        sid = row["Systematic ID"]
        draw_panel(
            ax,
            row,
            gallery.segments[sid],
            gallery.raw[sid],
            gallery.sites[sid],
            gallery.hand.get(sid, []),
            gallery.domains.get(sid),
            gallery.viability.get(sid),
            gallery.grna_dr.get(sid),
        )
        # every panel is labelled: a page is 6 stacked rows, and a reader who
        # screenshots one of them should not have to carry the axis meaning over
        # from a neighbour
        ax.set_xlabel("Residues(aa)")
        ax.set_ylabel("insertion DR")
    for ax in axes[len(page) :]:
        ax.set_visible(False)

    # the group heading moves up to the page: every panel now has its own title,
    # so a title on the first one would be overwritten by its gene's
    term = page["term"].iloc[0]
    fig = plt.gcf()
    fig.suptitle(
        f"{GENE_GROUP_LABELS[term]}  —  {start}-{end} / {total}",
        y=1 - SUPTITLE_STRIP / 2,
        va="center",
        fontweight="bold",
    )
    # hand fit_panels the legend column up front rather than letting it discover
    # the shortfall — see LEGEND_COLUMN_IN for what that discovery costs
    width, height = fig.get_size_inches()
    fig.set_size_inches(width + LEGEND_COLUMN_IN, height)
    fit_panels(rect=(0, 0, 1, 1 - SUPTITLE_STRIP))


@logger.catch(reraise=True)
def run(config: PlotConfig) -> int:
    """Render the gallery as one combined PDF plus one PDF per gene group."""
    apply_house_style()
    gallery = load_gallery(config)
    genes = gallery.genes

    # Pages never straddle a group: the heading names the group and the page's
    # position inside it, so a page holding two groups would have no one heading.
    groups: list[tuple[str, list[tuple[pd.DataFrame, int, int, int]]]] = []
    for term in (*TERM_ORDER, UNANNOTATED):
        subset = genes[genes["term"].eq(term)]
        total = len(subset)
        pages = [
            (subset.iloc[i : i + ROWS_PER_PAGE], i + 1, min(i + ROWS_PER_PAGE, total), total)
            for i in range(0, total, ROWS_PER_PAGE)
        ]
        groups.append((term, pages))

    def render(pages: list[tuple[pd.DataFrame, int, int, int]], stem: Path) -> None:
        def page_figures() -> Iterator[Figure]:
            for page, start, end, total in pages:
                draw_page(page, gallery, start, end, total)
                yield plt.gcf()
                plt.close()

        save_pdf_pages(stem, page_figures())

    all_pages = [page for _, pages in groups for page in pages]
    render(all_pages, config.output.with_suffix(""))
    # The count in the file name is what makes the split files recognisable in a
    # directory listing; the script, not Snakemake, is what knows it.
    split_dir = config.output.parent / config.output.stem
    # A count in the name means different data leaves the previous run's file
    # behind; remove only this script's own products, never the directory.
    for stale in (*split_dir.glob("*_genes.pdf"), *split_dir.glob("*_genes.review.png")):
        stale.unlink()
    for term, pages in groups:
        if pages:
            render(pages, split_dir / f"{GENE_GROUP_SLUGS[term]}_{pages[0][3]}_genes")

    counts = genes["term"].value_counts()
    logger.info(
        f"{len(genes):,} reported genes over {len(all_pages)} pages: "
        + ", ".join(f"{term} {count}" for term, count in counts.items())
    )
    logger.info(
        f"wrote {config.output.name}, its review PNG and "
        f"{sum(1 for _, pages in groups if pages)} group PDFs under {split_dir.name}/"
    )
    return 0


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Plot the domain-difference gallery (4a)")
    parser.add_argument("--segments", type=Path, required=True, help="Input domain_segments.tsv")
    parser.add_argument("--summary", type=Path, required=True, help="Input domain_summary.tsv")
    parser.add_argument(
        "--hand-regions", type=Path, required=True, help="Input domain_hand_regions.tsv"
    )
    parser.add_argument("--insertions", type=Path, required=True, help="Input insertions.parquet")
    parser.add_argument(
        "--annotation",
        type=Path,
        required=True,
        help="1c_annotation gene_annotation_reference parquet (panel titles)",
    )
    parser.add_argument(
        "--afdb-segments",
        type=Path,
        required=True,
        help="Upstream AFdb PAE/TED/ECOD segment table (resources/external copy)",
    )
    parser.add_argument("--output", type=Path, required=True, help="Output PDF path")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging"
    )
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render the gallery, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotConfig(
            segments=args.segments.resolve(),
            summary=args.summary.resolve(),
            hand_regions=args.hand_regions.resolve(),
            insertions=args.insertions.resolve(),
            annotation=args.annotation.resolve(),
            afdb_segments=args.afdb_segments.resolve(),
            output=args.output.resolve(),
        )
        config.output.parent.mkdir(parents=True, exist_ok=True)
        run(config)
    except (ValueError, FileNotFoundError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
