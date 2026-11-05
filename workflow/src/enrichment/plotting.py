"""
Enrichment Plotting
===================

Data preparation and rendering for the enrichment stage's two figure families:

- the per-variant GO dotplot (cluster x term, colour = significance, size = gene ratio)
- the cross-variant pair (a term x variant significance heatmap, and the count of
  significant terms per cluster for each variant)

Preparation functions are pure (frame in, frame out) and unit-tested; the draw
functions own the matplotlib axes. cnsplots supplies the house style and the
saving path (``figures.apply_house_style`` + ``figures.save_dual``), while the
bubble and heatmap marks themselves are raw Axes calls — cnsplots 0.7.0 has no
bubble chart, and its heatmap wrapper needs an AnnData plus a clustering tree
this figure does not have.

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
from pathlib import Path

# Local project
import cnsplots as cns
from loguru import logger
from matplotlib.axes import Axes

# Data processing / plotting
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from enrichment.schema import EnrichmentCol
from figures import DENSITY_CMAP, FURNITURE_COLOR, apply_house_style, house_colors, save_dual

# =============================================================================
# CONSTANTS
# =============================================================================
# Below this, -log10(FDR) saturates instead of returning inf for a p of exactly 0.
FDR_FLOOR = 1e-50
# Marker area in points^2, mapped linearly from the gene ratio.
MARKER_MIN_AREA = 12.0
MARKER_AREA_PER_RATIO = 260.0
# Row pitch of one term in the dotplot, in layout pixels (1 px = 1 pt).
DOTPLOT_ROW_PX = 11
DOTPLOT_WIDTH_PX = 320
DOTPLOT_MARGIN_PX = 96
DOTPLOT_COLUMNS = [
    EnrichmentCol.CLUSTER,
    EnrichmentCol.TERM_ID,
    EnrichmentCol.P_FDR,
    EnrichmentCol.GENE_RATIO,
]


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
def plot_dotplot(
    enrichment: pd.DataFrame,
    terms: pd.DataFrame,
    output_stem: Path,
    title: str,
    fdr_floor: float = FDR_FLOOR,
) -> None:
    """Render the cluster x term dotplot to ``output_stem``.pdf and .review.png."""
    matrix = dotplot_matrix(enrichment, terms, fdr_floor)
    n_rows = len(terms)
    marker_area = MARKER_MIN_AREA + MARKER_AREA_PER_RATIO * matrix[EnrichmentCol.GENE_RATIO]

    apply_house_style()
    cns.figure(width=DOTPLOT_WIDTH_PX, height=DOTPLOT_MARGIN_PX + DOTPLOT_ROW_PX * n_rows)
    figure = plt.gcf()
    ax = figure.subplots()

    scatter = ax.scatter(
        matrix[EnrichmentCol.CLUSTER].astype(int),
        matrix["row"],
        s=marker_area,
        c=matrix["neg_log10_fdr"],
        cmap=DENSITY_CMAP,
        linewidths=0.3,
        edgecolors=FURNITURE_COLOR,
    )

    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(terms[EnrichmentCol.TERM])
    ax.set_xticks(sorted(matrix[EnrichmentCol.CLUSTER].astype(int).unique()))
    ax.set_xlabel("Cluster  (dot size = gene ratio)")
    ax.set_ylabel("")
    ax.set_title(title, loc="left")
    ax.invert_yaxis()

    figure.colorbar(scatter, ax=ax, label=r"$-\log_{10}$(FDR)", pad=0.02, fraction=0.04)
    save_dual(output_stem)
    plt.close(figure)
    logger.success(f"Wrote dotplot -> {output_stem}.pdf")


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
    labels = [term_names.get(term_id, term_id) for term_id in matrix.index]
    n_rows = len(matrix)

    apply_house_style()
    cns.figure(width=DOTPLOT_WIDTH_PX, height=DOTPLOT_MARGIN_PX + DOTPLOT_ROW_PX * n_rows)
    figure = plt.gcf()
    ax = figure.subplots()

    image = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap=DENSITY_CMAP)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(labels)
    ax.set_title("Strongest enrichment per term", loc="left")
    figure.colorbar(image, ax=ax, label=r"$-\log_{10}$(FDR)", pad=0.02, fraction=0.04)

    save_dual(output_stem)
    plt.close(figure)
    logger.success(f"Wrote variant heatmap -> {output_stem}.pdf")


def plot_variant_counts(counts: pd.DataFrame, output_stem: Path, fdr_threshold: float) -> None:
    """Render grouped bars of significant-term counts per cluster, one bar per variant."""
    apply_house_style()
    cns.figure(width=DOTPLOT_WIDTH_PX, height=220)
    figure = plt.gcf()
    ax = figure.subplots()

    labels = list(counts.columns)
    colors = house_colors(range(len(labels)))
    clusters = np.arange(len(counts.index))
    width = 0.8 / len(labels)

    _draw_variant_bars(ax, counts, labels, colors, clusters, width)

    ax.set_xticks(clusters)
    ax.set_xticklabels(counts.index)
    ax.set_xlabel("Cluster")
    ax.set_ylabel(f"Significant terms (FDR < {fdr_threshold})")
    ax.legend(title="Variant", frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))

    save_dual(output_stem)
    plt.close(figure)
    logger.success(f"Wrote variant counts -> {output_stem}.pdf")


def _draw_variant_bars(
    ax: Axes,
    counts: pd.DataFrame,
    labels: list[str],
    colors: list[str],
    clusters: np.ndarray,
    width: float,
) -> None:
    """Draw one bar series per variant around each cluster position."""
    for offset, (label, color) in enumerate(zip(labels, colors, strict=True)):
        positions = clusters + (offset - (len(labels) - 1) / 2) * width
        ax.bar(positions, counts[label].to_numpy(), width=width, color=color, label=label)
