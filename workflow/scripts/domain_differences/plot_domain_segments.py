#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
- the grey band is the merge threshold ``lambda = z * sigma_gene + sigma_add``,
  centred on the gene's own level: a step that fits inside it is noise to this
  rule, so a segment only means something if its bar sits outside the band.
- green shading marks the hand-curated non-essential region, where there is one.

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

Output
------
- domain_segments.pdf — the gallery
- domain_segments.review.png — page 1, for a quick look

Usage
-----
    python workflow/scripts/domain_differences/plot_domain_segments.py \\
        --segments results/4a_domain_differences/HD_DIT_HAP/domain_segments.tsv \\
        --summary results/4a_domain_differences/HD_DIT_HAP/domain_summary.tsv \\
        --hand-regions results/4a_domain_differences/HD_DIT_HAP/domain_hand_regions.tsv \\
        --insertions results/4a_domain_differences/HD_DIT_HAP/_work/insertions.parquet \\
        --annotation results/1c_annotation/2026-06-01/2026-08-11/gene_annotation_reference.protein.parquet \\
        --output results/4a_domain_differences/HD_DIT_HAP/domain_segments.pdf

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-10-14
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Third-party Imports
from loguru import logger

# 3. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import cnsplots as cns  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from domain_differences.core import (  # noqa: E402
    TERM_ORDER,
    UNANNOTATED,
    collapse_sites,
    hand_region_map,
    read_gene_context,
)
from figures import (  # noqa: E402
    FURNITURE_COLOR,
    PanelShape,
    apply_house_style,
    fit_panels,
    grid_axes,
)
from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
ROWS_PER_PAGE = 6
GENE_GROUP_LABELS = {**{T: f"{T} terminus" for T in TERM_ORDER}, UNANNOTATED: "no hand annotation"}
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
    output: Path


@dataclass(kw_only=True, slots=True, frozen=True)
class Gallery:
    """Everything one render needs, grouped and indexed for lookup by gene."""

    genes: pd.DataFrame
    segments: dict
    sites: dict
    raw: dict
    hand: dict
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
    raw = {sid: frame for sid, frame in insertions.groupby("Systematic ID", sort=False)}
    sites = {sid: collapse_sites(frame) for sid, frame in raw.items()}
    viability, grna_dr = read_gene_context(config.annotation)

    # The term a gene is grouped under. No gene is annotated in two sheets, but if
    # that ever changed the first sheet in TERM_ORDER wins.
    term = (hand_table.dropna(subset=["sid"])
            .drop_duplicates("sid").set_index("sid")["term"].to_dict())
    rank = {t: i for i, t in enumerate((*TERM_ORDER, UNANNOTATED))}
    summary = pd.read_csv(config.summary, sep="\t")
    summary["term"] = [term.get(sid, UNANNOTATED) for sid in summary["Systematic ID"]]
    summary["_rank"] = summary["term"].map(rank).astype(int)
    genes = (summary.sort_values(["_rank", "max_boundary_contrast"], ascending=[True, False])
             .drop(columns="_rank").reset_index(drop=True))

    return Gallery(
        genes=genes,
        segments={sid: group for sid, group in segments.groupby("Systematic ID", sort=False)},
        sites=sites, raw=raw, hand=hand,
        viability=viability, grna_dr=grna_dr,
    )


# =============================================================================
# CORE LOGIC — one panel
# =============================================================================
def draw_insertions(ax, raw: pd.DataFrame, sites: pd.DataFrame) -> None:
    """Raw Forward/Reverse points of one gene plus their per-site mean.

    The typical Forward-vs-Reverse gap is only 0.045 DR, a couple of pixels at this
    page scale, so a pair drawn as two dots looks like one point. The thin vertical
    tick joining them makes every site read as a pair, and its length IS the
    disagreement. The mean is drawn as an OPEN diamond so it never hides the pair.
    """
    forward = raw["Insertion_direction"].eq("Forward").to_numpy()
    residues = raw["Residue_affected"].astype(int).to_numpy()
    ax.vlines(sites["Residue_affected"].astype(int),
              sites["DR_rev"].fillna(sites["DR"]), sites["DR_fwd"].fillna(sites["DR"]),
              color="0.6", lw=0.7, zorder=3)
    ax.scatter(residues[forward], raw["DR"].to_numpy()[forward],
               s=9, color=FWD_COLOR, lw=0, zorder=4)
    ax.scatter(residues[~forward], raw["DR"].to_numpy()[~forward],
               s=9, color=REV_COLOR, lw=0, zorder=4)
    ax.scatter(sites["Residue_affected"].astype(int), sites["DR"], s=16, marker="D",
               facecolors="none", edgecolors=SITE_MEAN_COLOR, linewidths=0.8, zorder=5)


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
    call = viability if viability else "—"
    # `grna_dr == grna_dr` is False only for NaN, which read_gene_context already
    # filters out — kept because a NaN that printed as "gRNA DR=+nan" would look
    # like a value rather than a missing one
    dr = f"gRNA DR={grna_dr:+.2f}" if grna_dr is not None and grna_dr == grna_dr else "gRNA DR=—"
    return f"{name}  —  {call}\n{dr}"


def draw_panel(ax, row: pd.Series, segments: pd.DataFrame, raw: pd.DataFrame,
               sites: pd.DataFrame, hand: list[tuple[int, int]],
               viability: str | None = None, grna_dr: float | None = None) -> None:
    """One gene: heading, threshold band, hand shading, the two channels, segment levels, caption."""
    # Sized so the gallery's widest caption still fits its panel: at 6.5 pt the
    # longest line (a systematic-ID name plus the three cap fields) runs 149 pt
    # into a 150 pt panel and crosses the spine, and 6.0 leaves 9 pt of margin.
    caption_size = plt.rcParams["font.size"] - 2.0
    ax.set_title(panel_title(str(row["name"]), viability, grna_dr),
                 fontsize=plt.rcParams["font.size"] - 1.0)
    ax.axhline(0, color=LAMBDA_COLOR, lw=0.6, zorder=0)

    lam = float(row["lambda"])
    level = float(sites["DR"].mean())
    if lam > 0:
        ax.axhspan(level - lam / 2, level + lam / 2, color=LAMBDA_COLOR, alpha=0.13, lw=0, zorder=1)
        for edge in (level - lam / 2, level + lam / 2):
            ax.axhline(edge, color=LAMBDA_COLOR, lw=0.8, ls=":", zorder=2)
        ax.text(0.98, level + lam / 2, f"λ={lam:.2f}", transform=ax.get_yaxis_transform(),
                ha="right", va="bottom", color=SITE_MEAN_COLOR, fontsize=caption_size)

    for start, end in hand:
        ax.axvspan(start, end, color=HAND_COLOR, alpha=0.18, lw=0, zorder=1)

    draw_insertions(ax, raw, sites)
    for _, segment in segments.iterrows():
        ax.hlines(float(segment["mean_DR"]), int(segment["start_residue"]),
                  int(segment["end_residue"]), color=SEGMENT_COLOR, lw=1.6, zorder=6)
    for _, segment in segments.iloc[1:].iterrows():
        ax.axvline(int(segment["start_residue"]), color="0.55", lw=0.8, ls="--", zorder=5)

    # x range covers the hand region as well as the insertions: an annotated terminus
    # no insertion reached is exactly the case a reader must be able to see
    bounds = [float(sites["Residue_affected"].min()), float(sites["Residue_affected"].max())]
    for start, end in hand:
        bounds += [float(start), float(end)]
    pad = max(0.02 * max(bounds), 1.0)
    ax.set_xlim(min(bounds) - pad, max(bounds) + pad)
    # see CAPTION_BAND: the data window ends up occupying the top 1 - CAPTION_BAND
    # of the axes, whatever the panel's absolute size turns out to be
    low, high = Y_RANGE
    ax.set_ylim(low - CAPTION_BAND / (1 - CAPTION_BAND) * (high - low), high)

    ranges = "|".join(f"{int(s['start_residue'])}-{int(s['end_residue'])}"
                      for _, s in segments.iloc[:SEGMENTS_SHOWN].iterrows())
    if len(segments) > SEGMENTS_SHOWN:
        ranges += f"|+{len(segments) - SEGMENTS_SHOWN}"
    ax.text(0.02, 0.04,
            f"{row['name']}  segs={int(row['n_segments'])}  n_sites={int(row['n_sites'])}"
            f"  σ={float(row['sigma_gene']):.2f}  λ={lam:.2f}\n"
            f"{ranges}  max|Δ|={float(row['max_boundary_contrast']):.2f}",
            transform=ax.transAxes, va="bottom", linespacing=1.3, fontsize=caption_size)

    # every row carries the key, in one column beside it: a reader who scrolled to
    # a single gene should not have to find the page's legend to know which dots
    # are Forward. Outside the axes, so it never sits on the data.
    ax.legend(handles=legend_handles(), loc="center left", bbox_to_anchor=(1.02, 0.5),
              ncol=1, frameon=False, fontsize=caption_size, handletextpad=0.6,
              borderaxespad=0.0, labelspacing=0.8)


def legend_handles() -> list[Line2D]:
    """One key per visual channel, with the site mean open like it is in the panels."""
    return [
        Line2D([], [], marker="o", ls="", color=FWD_COLOR, label="insertion Forward"),
        Line2D([], [], marker="o", ls="", color=REV_COLOR, label="insertion Reverse"),
        Line2D([], [], marker="D", ls="", markerfacecolor="none",
               markeredgecolor=SITE_MEAN_COLOR, label="site mean (segmented on)"),
        Line2D([], [], color=SEGMENT_COLOR, lw=1.6, label="segment mean DR"),
        Line2D([], [], color=LAMBDA_COLOR, lw=6, alpha=0.35, label="λ = z·σ + c band"),
        Line2D([], [], color=HAND_COLOR, lw=6, alpha=0.35, label="hand-curated non-essential"),
    ]


# =============================================================================
# CORE LOGIC — the gallery
# =============================================================================
def draw_page(page: pd.DataFrame, gallery: Gallery) -> None:
    """Draw one page of up to ROWS_PER_PAGE genes into a freshly created figure."""
    # empty labels: a row is a named gene, not a lettered panel, and its own title
    # already carries the identity
    axes = grid_axes(ROWS_PER_PAGE, 1, labels=[], shape=PanelShape.WIDE)
    for ax, (_, row) in zip(axes, page.iterrows()):
        sid = row["Systematic ID"]
        draw_panel(ax, row, gallery.segments[sid], gallery.raw[sid], gallery.sites[sid],
                   gallery.hand.get(sid, []),
                   gallery.viability.get(sid), gallery.grna_dr.get(sid))
        # every panel is labelled: a page is 6 stacked rows, and a reader who
        # screenshots one of them should not have to carry the axis meaning over
        # from a neighbour
        ax.set_xlabel("Residues(aa)")
        ax.set_ylabel("insertion DR")
    for ax in axes[len(page):]:
        ax.set_visible(False)

    # the group heading moves up to the page: every panel now has its own title,
    # so a title on the first one would be overwritten by its gene's
    term = page["term"].iloc[0]
    fig = plt.gcf()
    fig.suptitle(f"{GENE_GROUP_LABELS[term]}  —  {len(page)} genes",
                 y=1 - SUPTITLE_STRIP / 2, va="center",
                 fontsize=plt.rcParams["font.size"] + 2.0, fontweight="bold")
    # hand fit_panels the legend column up front rather than letting it discover
    # the shortfall — see LEGEND_COLUMN_IN for what that discovery costs
    width, height = fig.get_size_inches()
    fig.set_size_inches(width + LEGEND_COLUMN_IN, height)
    fit_panels(rect=(0, 0, 1, 1 - SUPTITLE_STRIP))


@logger.catch(reraise=True)
def run(config: PlotConfig) -> int:
    """Render every page of the gallery into one PDF, plus page 1 as a preview PNG."""
    apply_house_style()
    gallery = load_gallery(config)
    genes = gallery.genes
    pages = [genes.iloc[i:i + ROWS_PER_PAGE] for i in range(0, len(genes), ROWS_PER_PAGE)]
    review = config.output.with_suffix(".review.png")

    with PdfPages(config.output) as pdf:
        for number, page in enumerate(pages, start=1):
            draw_page(page, gallery)
            figure = plt.gcf()
            pdf.savefig(figure, dpi=300, bbox_inches="tight")
            if number == 1:
                # cns.savefig writes only the current figure, so the page-1 preview
                # is saved here rather than through save_dual(), whose PDF half
                # would clobber the gallery PdfPages is still writing
                cns.savefig(review)
            plt.close(figure)

    counts = genes["term"].value_counts()
    logger.info(f"{len(genes):,} reported genes over {len(pages)} pages: "
                + ", ".join(f"{term} {count}" for term, count in counts.items()))
    logger.info(f"wrote {config.output.name} and {review.name}")
    return 0


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Plot the domain-difference gallery (4a)")
    parser.add_argument("--segments", type=Path, required=True, help="Input domain_segments.tsv")
    parser.add_argument("--summary", type=Path, required=True, help="Input domain_summary.tsv")
    parser.add_argument("--hand-regions", type=Path, required=True, help="Input domain_hand_regions.tsv")
    parser.add_argument("--insertions", type=Path, required=True, help="Input insertions.parquet")
    parser.add_argument("--annotation", type=Path, required=True,
                        help="1c_annotation gene_annotation_reference parquet (panel titles)")
    parser.add_argument("--output", type=Path, required=True, help="Output PDF path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render the gallery, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotConfig(
            segments=args.segments.resolve(), summary=args.summary.resolve(),
            hand_regions=args.hand_regions.resolve(), insertions=args.insertions.resolve(),
            annotation=args.annotation.resolve(), output=args.output.resolve(),
        )
        config.output.parent.mkdir(parents=True, exist_ok=True)
        run(config)
    except (ValueError, FileNotFoundError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
