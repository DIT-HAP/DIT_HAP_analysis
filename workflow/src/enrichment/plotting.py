"""
Enrichment Plotting
===================

Data preparation and rendering for the enrichment stage's two figure families:

- the per-variant GO dotplot (cluster x term, colour = significance, size = gene ratio)
- the cross-variant pair (a term x variant significance heatmap, and the count of
  significant terms per cluster for each variant)

Preparation functions are pure (frame in, frame out) and unit-tested; the draw
functions own the matplotlib axes. The dotplot matrix is drawn by ``cns.dotplot``
(the house enrichment-by-size-and-colour API). The cross-variant heatmap and
counts stay raw Axes calls: ``cns.heatmapplot`` accepts a plain DataFrame, but it
puts its colour legend outside the axes where it collides with the title at this
page budget, and ``cns.barplot`` draws its key inside the axes, over the bars.
An ``imshow`` + colourbar and a grouped ``bar`` keep both keys inside the space
the canvas reserves for them.

Every canvas is budgeted in points against the house 100 pt panel contract:
``cns.figure``'s "pixel" is 1/72 inch, so its numbers are points, and the page
grows around a fixed-size matrix rather than the matrix shrinking to fit a
guessed page. Panel boxes come from the house shapes and tick sizes from the
resolved cnsplots style, so no font or marker size is set here.

Input
-----
- ``go_enrichment_full_filtered.tsv`` (one per variant) from finalize_enrichment

Output
------
- PDF + ``.review.png`` pairs, written by ``figures.save_dual``
"""

# =============================================================================
# IMPORTS
# =============================================================================
# Standard library
from collections.abc import Iterable
from pathlib import Path

# Local project
import cnsplots as cns
from loguru import logger
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties

# Data processing / plotting
import matplotlib.pyplot as plt
from matplotlib.textpath import TextToPath
import numpy as np
import pandas as pd

from enrichment.schema import EnrichmentCol
from figures import (
    DENSITY_CMAP,
    PANEL_SIZES_PX,
    PanelShape,
    apply_house_style,
    house_colors,
    save_dual,
)

# =============================================================================
# CONSTANTS
# =============================================================================
# Below this, -log10(FDR) saturates instead of returning inf for a p of exactly 0.
FDR_FLOOR = 1e-50
# Every canvas measurement is in points, the unit cnsplots sizes pages in
# (cns.figure divides by 72 to get inches) and the unit the 100 pt panel contract
# is written in. Panel boxes come from the house shapes; tick sizes are read back
# from the resolved cnsplots style at draw time rather than restated here.
MATRIX_PT = PANEL_SIZES_PX[PanelShape.SQUARE][0]  # the dot matrix or heatmap block
ROW_GUTTER_PT = 3  # air between one categorical row's tick label and the next
# cns.dotplot measures its legend area in millimetres, not points; 20 mm is the
# API default and the width this canvas reserves for it.
LEGEND_MM = 20
LABEL_PAD_PT = 4  # gutter left of the longest row label
TICK_PAD_PT = 2
RIGHT_PAD_PT = 8  # the size legend overruns its own axes by ~6 pt
TITLE_PT = 18  # title band above the matrix
FOOTER_PT = 40  # x tick labels, their 15 pt labelpad and the xlabel
HEATMAP_TICK_ANGLE_DEG = 45  # variant names do not fit on a 12 pt column unrotated
# Hand-placed colorbar column: gap, bar, then its tick labels and axis label.
COLORBAR_GAP_PT = 6
COLORBAR_BAR_PT = 8
COLORBAR_TEXT_PT = 26
# Grouped bars need a wide panel: nine cluster groups of one bar per variant.
COUNTS_LEFT_PT = 40  # y tick labels plus the rotated ylabel
COUNTS_MATRIX_WIDTH_PT, COUNTS_MIN_HEIGHT_PT = PANEL_SIZES_PX[PanelShape.WIDE]
COUNTS_LEGEND_PAD_PT = 30  # key marker, its handletext pad, and the key title
COUNTS_FOOTER_PT = 30  # x tick labels plus the xlabel
COUNTS_TOP_PT = 8
PT_PER_MM = 72 / 25.4
# Column names double as the cns.dotplot legend titles, so they are written for
# display rather than as the tidy frame's identifiers.
DOTPLOT_COLOR_COLUMN = "-log10(FDR)"
DOTPLOT_SIZE_COLUMN = "Gene ratio"
DOTPLOT_COLUMNS = [
    EnrichmentCol.CLUSTER,
    EnrichmentCol.TERM_ID,
    EnrichmentCol.P_FDR,
    EnrichmentCol.GENE_RATIO,
]


# =============================================================================
# SHARED MEASUREMENT
# =============================================================================
def _tick_font_pt() -> float:
    """The tick size the cnsplots style resolved, which the canvas must budget."""
    resolved = FontProperties(size=plt.rcParams["xtick.labelsize"])
    return float(resolved.get_size_in_points())


def _row_pitch_pt() -> float:
    """Vertical pitch of one categorical row: its tick label plus the gutter."""
    return _tick_font_pt() + ROW_GUTTER_PT


# =============================================================================
# DOTPLOT — DATA
# =============================================================================
def dotplot_terms(enrichment: pd.DataFrame, top_terms_per_cluster: int) -> pd.DataFrame:
    """Return the union of each cluster's most significant terms, grouped by home cluster.

    A term's home cluster is the one where it reaches its smallest p-value. Ordering
    rows by home cluster keeps each block's dots in one column, so the panel reads as
    blocks (``cluster`` then ``namespace`` then significance) instead of a scatter.
    """
    ranked = enrichment.sort_values([EnrichmentCol.CLUSTER, EnrichmentCol.P_FDR])
    top = ranked.groupby(EnrichmentCol.CLUSTER, sort=True).head(top_terms_per_cluster)

    columns = [EnrichmentCol.NAMESPACE, EnrichmentCol.TERM_ID, EnrichmentCol.TERM]
    terms = top[columns].drop_duplicates(EnrichmentCol.TERM_ID).copy()

    strongest = top.loc[top.groupby(EnrichmentCol.TERM_ID)[EnrichmentCol.P_FDR].idxmin()]
    strongest = strongest.set_index(EnrichmentCol.TERM_ID)
    terms["home_cluster"] = terms[EnrichmentCol.TERM_ID].map(strongest[EnrichmentCol.CLUSTER])
    terms["best_p"] = terms[EnrichmentCol.TERM_ID].map(strongest[EnrichmentCol.P_FDR])
    return terms.sort_values(
        ["home_cluster", EnrichmentCol.NAMESPACE, "best_p", EnrichmentCol.TERM_ID]
    )


def dotplot_matrix(
    enrichment: pd.DataFrame, terms: pd.DataFrame, fdr_floor: float = FDR_FLOOR
) -> pd.DataFrame:
    """Select the rows for the chosen terms and add the plotting columns."""
    is_selected = enrichment[EnrichmentCol.TERM_ID].isin(terms[EnrichmentCol.TERM_ID])
    selected = enrichment.loc[is_selected, DOTPLOT_COLUMNS].copy()
    selected["neg_log10_fdr"] = -np.log10(selected[EnrichmentCol.P_FDR].clip(lower=fdr_floor))
    selected["row"] = selected[EnrichmentCol.TERM_ID].map(
        {term_id: row for row, term_id in enumerate(terms[EnrichmentCol.TERM_ID])}
    )
    return selected


# =============================================================================
# DOTPLOT — RENDER
# =============================================================================
def dotplot_plot_frame(
    enrichment: pd.DataFrame, terms: pd.DataFrame, fdr_floor: float = FDR_FLOOR
) -> pd.DataFrame:
    """Shape the selected cells into the x / y-category / colour / size frame cns.dotplot reads."""
    matrix = dotplot_matrix(enrichment, terms, fdr_floor)
    clusters = matrix[EnrichmentCol.CLUSTER].astype(int)
    names = terms.set_index(EnrichmentCol.TERM_ID)[EnrichmentCol.TERM]
    return pd.DataFrame(
        {
            # An ordered categorical keeps cluster 10 after cluster 9; the plotter
            # pivots lexicographically on a plain string column.
            EnrichmentCol.CLUSTER: pd.Categorical(
                clusters.astype(str),
                categories=[str(cluster) for cluster in sorted(clusters.unique())],
                ordered=True,
            ),
            "term": pd.Categorical(
                matrix[EnrichmentCol.TERM_ID].map(names),
                categories=list(terms[EnrichmentCol.TERM]),
                ordered=True,
            ),
            DOTPLOT_COLOR_COLUMN: matrix["neg_log10_fdr"],
            DOTPLOT_SIZE_COLUMN: matrix[EnrichmentCol.GENE_RATIO],
        }
    )


def _installed_house_font() -> str:
    """The first face of the house stack this machine actually has.

    Measuring against the whole fallback chain makes matplotlib log one
    "font family not found" warning per absent face, which is noise: only the
    first installed face is ever used for measurement.
    """
    installed = {font.name for font in font_manager.fontManager.ttflist}
    return next(name for name in plt.rcParams["font.sans-serif"] if name in installed)


def label_width_pt(names: Iterable[str]) -> float:
    """Width the longest label needs at the house tick size, in points."""
    font = FontProperties(family=_installed_house_font(), size=_tick_font_pt())
    engine = TextToPath()
    widths = [engine.get_text_width_height_descent(str(name), font, False)[0] for name in names]
    return max(widths)


def rotated_label_height_pt(names: Iterable[str], angle_deg: float) -> float:
    """Height a band of tick labels rotated by ``angle_deg`` needs, in points."""
    angle = np.radians(angle_deg)
    width_pt = label_width_pt(names)
    return float(width_pt * np.sin(angle) + _tick_font_pt() * np.cos(angle))


def plot_dotplot(
    enrichment: pd.DataFrame,
    terms: pd.DataFrame,
    output_stem: Path,
    title: str,
    fdr_floor: float = FDR_FLOOR,
) -> None:
    """Render the cluster x term dotplot to ``output_stem``.pdf and .review.png."""
    apply_house_style()

    if terms.empty:
        _plot_empty_dotplot(output_stem, title)
        return

    frame = dotplot_plot_frame(enrichment, terms, fdr_floor)
    legend_pt = LEGEND_MM * PT_PER_MM
    label_col_pt = LABEL_PAD_PT + TICK_PAD_PT + label_width_pt(terms[EnrichmentCol.TERM])
    # The canvas grows around a fixed 100 pt matrix: the plotter reserves the legend
    # area inside the host axes, so the matrix keeps its size whatever the labels cost.
    width = round(label_col_pt + MATRIX_PT + legend_pt + RIGHT_PAD_PT)
    height = round(TITLE_PT + _row_pitch_pt() * len(terms) + FOOTER_PT)

    cns.figure(width=width, height=height, unit="pt")
    figure = plt.gcf()
    # add_axes, not subplots: pycomplexheatmap's colorbar calls subplots_adjust(bottom=0),
    # which would wipe out the footer this canvas budgets for.
    host = figure.add_axes(
        (
            label_col_pt / width,
            FOOTER_PT / height,
            (MATRIX_PT + legend_pt) / width,
            1 - (TITLE_PT + FOOTER_PT) / height,
        )
    )

    # value= collapses the pivot to a single numeric column: without it the
    # plotter aggregates every numeric column and pycomplexheatmap's stack()
    # yields four. The *_na=NaN fills keep the matrix sparse, so an untested
    # cluster-by-term cell is blank instead of a dot at the zero end of the scale.
    plotter = cns.dotplot(
        data=frame,
        x=EnrichmentCol.CLUSTER,
        y="term",
        color=DOTPLOT_COLOR_COLUMN,
        size=DOTPLOT_SIZE_COLUMN,
        value=DOTPLOT_COLOR_COLUMN,
        value_na=np.nan,
        s_na=np.nan,
        c_na=np.nan,
        ax=host,
        xticklabels_rotation=0,
        # extend="neither" squares off the colorbar ends the plotter rounds by default;
        # %.0f drops the three decimals its auto ticks carry at this data range.
        cmap_legend_kws={"extend": "neither", "format": "%.0f"},
        xlabel="Cluster  (dot size = gene ratio)",
        ylabel="",
        cmap=DENSITY_CMAP,
    )
    plotter.hm_ax.set_title(title, loc="left")
    # cnsplots right-aligns x tick labels for its rotated default; at rotation 0
    # that leaves every label half an advance to the left of its tick.
    plt.setp(plotter.hm_ax.get_xticklabels(), ha="center")

    save_dual(output_stem)
    plt.close(figure)
    logger.success(f"Wrote dotplot -> {output_stem}.pdf")


def _plot_empty_dotplot(output_stem: Path, title: str) -> None:
    """Write a one-line note figure when a variant has no significant term to plot."""
    width = LABEL_PAD_PT + MATRIX_PT + RIGHT_PAD_PT
    height = round(TITLE_PT + _row_pitch_pt() + FOOTER_PT)
    cns.figure(width=width, height=height, unit="pt")
    figure = plt.gcf()
    ax = figure.add_axes((0, 0, 1, 1))
    ax.axis("off")
    ax.set_title(title, loc="left")
    ax.text(0.5, 0.5, "no significant terms", ha="center", va="center")

    save_dual(output_stem)
    plt.close(figure)
    logger.warning(f"No significant terms; wrote a placeholder dotplot -> {output_stem}.pdf")


# =============================================================================
# CROSS-VARIANT — DATA
# =============================================================================
def variant_heatmap_matrix(
    per_variant: dict[str, pd.DataFrame], top_terms: int, fdr_floor: float = FDR_FLOOR
) -> pd.DataFrame:
    """Best -log10(FDR) per term and variant, over the union of each variant's top terms.

    Missing cells stay NaN: a term that is not significant in a variant is a real
    absence, not a zero, so the cell is left blank rather than painted as a miss.
    Rows are ordered by namespace, then by the strongest value reached in any
    variant, so the heatmap reads as blocks instead of an alphabetical run.
    """
    best_p = {
        label: frame.groupby(EnrichmentCol.TERM_ID)[EnrichmentCol.P_FDR].min()
        for label, frame in per_variant.items()
    }
    candidates: set[str] = {
        str(term_id)
        for per_term in best_p.values()
        for term_id in per_term.nsmallest(top_terms).index
    }
    selected_terms = pd.Index(sorted(candidates))

    matrix = pd.DataFrame(
        {label: per_term.reindex(selected_terms) for label, per_term in best_p.items()},
        index=selected_terms,
    )
    clipped = matrix.clip(lower=fdr_floor)
    scores = pd.DataFrame(
        -np.log10(clipped.to_numpy(dtype=float)), index=clipped.index, columns=clipped.columns
    )
    return scores.loc[_heatmap_row_order(scores, per_variant)]


def _heatmap_row_order(scores: pd.DataFrame, per_variant: dict[str, pd.DataFrame]) -> pd.Index:
    """Order terms by namespace, then by their strongest value across variants."""
    namespaces = pd.concat(
        [frame[[EnrichmentCol.TERM_ID, EnrichmentCol.NAMESPACE]] for frame in per_variant.values()]
    ).drop_duplicates(EnrichmentCol.TERM_ID)
    namespaces = namespaces.set_index(EnrichmentCol.TERM_ID)[EnrichmentCol.NAMESPACE]

    ranked = pd.DataFrame(
        {
            EnrichmentCol.NAMESPACE: namespaces.reindex(scores.index),
            "best": scores.max(axis=1),
        },
        index=scores.index,
    )
    return ranked.sort_values([EnrichmentCol.NAMESPACE, "best"], ascending=[True, False]).index


def significant_term_counts(
    per_variant: dict[str, pd.DataFrame], fdr_threshold: float
) -> pd.DataFrame:
    """Count FDR-significant terms per cluster for every variant, zero-filled."""
    counts: dict[str, pd.Series] = {}
    for label, frame in per_variant.items():
        is_significant = frame[EnrichmentCol.P_FDR] < fdr_threshold
        counts[label] = frame.loc[is_significant].groupby(EnrichmentCol.CLUSTER).size()
    return pd.DataFrame(counts).fillna(0).astype(int).sort_index()


def term_display_names(per_variant: dict[str, pd.DataFrame]) -> pd.Series:
    """Map every term id seen in any variant to its display name."""
    named = pd.concat(
        [frame[[EnrichmentCol.TERM_ID, EnrichmentCol.TERM]] for frame in per_variant.values()]
    ).drop_duplicates(EnrichmentCol.TERM_ID)
    return named.set_index(EnrichmentCol.TERM_ID)[EnrichmentCol.TERM]


# =============================================================================
# CROSS-VARIANT — RENDER
# =============================================================================
def plot_variant_heatmap(matrix: pd.DataFrame, term_names: pd.Series, output_stem: Path) -> None:
    """Render the term x variant significance heatmap to ``output_stem``.pdf and .review.png."""
    apply_house_style()

    labels = term_names.loc[matrix.index].tolist()
    footer_pt = rotated_label_height_pt(matrix.columns, HEATMAP_TICK_ANGLE_DEG) + TICK_PAD_PT
    label_col_pt = LABEL_PAD_PT + TICK_PAD_PT + label_width_pt(labels)
    width = round(label_col_pt + MATRIX_PT + COLORBAR_GAP_PT + COLORBAR_BAR_PT + COLORBAR_TEXT_PT)
    height = round(TITLE_PT + _row_pitch_pt() * len(matrix) + footer_pt)

    cns.figure(width=width, height=height, unit="pt")
    figure = plt.gcf()
    matrix_left = label_col_pt / width
    matrix_bottom = footer_pt / height
    matrix_height_frac = 1 - (TITLE_PT + footer_pt) / height
    ax = figure.add_axes((matrix_left, matrix_bottom, MATRIX_PT / width, matrix_height_frac))

    image = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap=DENSITY_CMAP)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=HEATMAP_TICK_ANGLE_DEG, ha="right")
    ax.set_yticks(range(len(matrix)))
    ax.set_yticklabels(labels)
    ax.set_title("Strongest enrichment per term", loc="left")

    cbar_left = matrix_left + (MATRIX_PT + COLORBAR_GAP_PT) / width
    cax = figure.add_axes((cbar_left, matrix_bottom, COLORBAR_BAR_PT / width, matrix_height_frac))
    figure.colorbar(image, cax=cax, label=r"$-\log_{10}$(FDR)", extend="neither")

    save_dual(output_stem)
    plt.close(figure)
    logger.success(f"Wrote variant heatmap -> {output_stem}.pdf")


def plot_variant_counts(counts: pd.DataFrame, output_stem: Path, fdr_threshold: float) -> None:
    """Render grouped bars of significant-term counts per cluster, one bar per variant."""
    apply_house_style()

    labels = list(counts.columns)
    ylabel = f"Significant terms (FDR < {fdr_threshold})"
    # The key sits outside the axes and the ylabel is rotated, so both can outgrow
    # the panel when a variant name is long: measure them instead of guessing.
    legend_pt = label_width_pt(labels) + COUNTS_LEGEND_PAD_PT
    axes_height_pt = max(COUNTS_MIN_HEIGHT_PT, label_width_pt([ylabel]) + TICK_PAD_PT)
    width = round(COUNTS_LEFT_PT + COUNTS_MATRIX_WIDTH_PT + legend_pt)
    height = round(COUNTS_TOP_PT + axes_height_pt + COUNTS_FOOTER_PT)

    cns.figure(width=width, height=height, unit="pt")
    figure = plt.gcf()
    ax = figure.add_axes(
        (
            COUNTS_LEFT_PT / width,
            COUNTS_FOOTER_PT / height,
            COUNTS_MATRIX_WIDTH_PT / width,
            axes_height_pt / height,
        )
    )

    colors = house_colors(range(len(labels)))
    clusters = np.arange(len(counts.index))
    width_per_bar = 0.8 / len(labels)
    for offset, (label, color) in enumerate(zip(labels, colors, strict=True)):
        positions = clusters + (offset - (len(labels) - 1) / 2) * width_per_bar
        ax.bar(positions, counts[label].to_numpy(), width=width_per_bar, color=color, label=label)

    ax.set_xticks(clusters)
    ax.set_xticklabels(counts.index)
    ax.set_xlabel("Cluster")
    ax.set_ylabel(ylabel)
    ax.legend(title="Variant", loc="upper left", bbox_to_anchor=(1.02, 1.0))

    save_dual(output_stem)
    plt.close(figure)
    logger.success(f"Wrote variant counts -> {output_stem}.pdf")
