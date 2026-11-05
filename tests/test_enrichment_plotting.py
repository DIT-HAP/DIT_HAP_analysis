"""Tests for the pure data preparation behind the enrichment figures."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.src.enrichment.plotting import (
    DOTPLOT_COLOR_COLUMN,
    DOTPLOT_SIZE_COLUMN,
    dotplot_matrix,
    dotplot_plot_frame,
    dotplot_terms,
    label_width_pt,
    rotated_label_height_pt,
    significant_term_counts,
    term_display_names,
    variant_heatmap_matrix,
)

COLUMNS = ["Cluster", "namespace", "term_id", "term", "p_fdr", "gene_ratio"]


def _enrichment(rows: list[list]) -> pd.DataFrame:
    """Build a filtered-GO frame with the columns the plotting helpers read."""
    return pd.DataFrame(rows, columns=COLUMNS)


def _two_cluster_frame() -> pd.DataFrame:
    return _enrichment(
        [
            [1, "BP", "GO:1", "process A", 1e-5, 0.10],
            [1, "BP", "GO:2", "process B", 1e-3, 0.20],
            [1, "CC", "GO:3", "component C", 1e-2, 0.30],
            [2, "BP", "GO:2", "process B", 1e-6, 0.25],
            [2, "CC", "GO:4", "component D", 1e-4, 0.15],
        ]
    )


def test_dotplot_terms_takes_the_union_of_per_cluster_tops():
    """Two terms per cluster -> the union of both clusters' top two."""
    terms = dotplot_terms(_two_cluster_frame(), top_terms_per_cluster=2)
    assert set(terms["term_id"]) == {"GO:1", "GO:2", "GO:4"}


def test_dotplot_terms_orders_by_home_cluster_then_namespace():
    """Rows group by the cluster where each term is strongest; namespace breaks ties."""
    terms = dotplot_terms(_two_cluster_frame(), top_terms_per_cluster=2)
    assert list(terms["term_id"]) == ["GO:1", "GO:2", "GO:4"]
    assert list(terms["home_cluster"]) == [1, 2, 2]


def test_dotplot_matrix_maps_rows_and_colours():
    """Every selected row carries its term row index and a -log10(FDR) value."""
    enrichment = _two_cluster_frame()
    terms = dotplot_terms(enrichment, top_terms_per_cluster=2)
    matrix = dotplot_matrix(enrichment, terms)

    assert len(matrix) == 4  # (cluster, term) pairs kept from the selected terms
    assert matrix["row"].max() == len(terms) - 1
    strongest = matrix.loc[matrix["term_id"] == "GO:2", "neg_log10_fdr"].max()
    assert strongest == pytest.approx(6.0)


def test_dotplot_plot_frame_is_sparse_and_keeps_term_order():
    """cns.dotplot gets one row per tested cell, with the term axis in home-cluster order."""
    enrichment = _two_cluster_frame()
    terms = dotplot_terms(enrichment, top_terms_per_cluster=2)
    frame = dotplot_plot_frame(enrichment, terms)

    assert list(frame.columns) == [
        "Cluster",
        "term",
        DOTPLOT_COLOR_COLUMN,
        DOTPLOT_SIZE_COLUMN,
    ]
    assert len(frame) == 4  # 3 terms x 2 clusters = 6 cells, only the 4 tested ones
    assert list(frame["term"].cat.categories) == [
        "process A",
        "process B",
        "component D",
    ]
    assert set(frame["Cluster"]) == {"1", "2"}


def test_dotplot_plot_frame_orders_clusters_numerically():
    """Cluster order is numeric, so cluster 10 follows cluster 9 instead of cluster 1."""
    enrichment = _enrichment(
        [
            [cluster, "BP", f"GO:{cluster}", f"process {cluster}", 1e-3, 0.1]
            for cluster in range(1, 11)
        ]
    )
    terms = dotplot_terms(enrichment, top_terms_per_cluster=1)
    frame = dotplot_plot_frame(enrichment, terms)

    assert list(frame["Cluster"].cat.categories) == [
        str(cluster) for cluster in range(1, 11)
    ]


def test_variant_heatmap_leaves_unseen_cells_empty():
    """A term absent from a variant stays NaN rather than becoming a zero."""
    first = _two_cluster_frame()
    second = _enrichment([[1, "BP", "GO:1", "process A", 1e-8, 0.10]])
    matrix = variant_heatmap_matrix({"first": first, "second": second}, top_terms=2)

    assert pd.isna(matrix.loc["GO:2", "second"])
    assert matrix.loc["GO:1", "second"] == pytest.approx(8.0)


def test_significant_term_counts_counts_per_cluster_and_variant():
    """Counts are per (variant, cluster) and zero-filled where a cluster has none."""
    counts = significant_term_counts(
        {
            "first": _two_cluster_frame(),
            "second": _enrichment([[1, "BP", "GO:1", "a", 0.5, 0.1]]),
        },
        fdr_threshold=1e-2,
    )
    assert counts.loc[1, "first"] == 2
    assert counts.loc[2, "first"] == 2
    assert counts.loc[1, "second"] == 0


def test_term_display_names_maps_every_seen_term():
    """The name map covers terms from both variants."""
    names = term_display_names(
        {"first": _two_cluster_frame(), "second": _two_cluster_frame()}
    )
    assert names["GO:4"] == "component D"


def test_label_width_pt_measures_the_longest_name():
    """The canvas reserves room for the widest label, so the longest name must win."""
    short = label_width_pt(pd.Series(["protein transport", "ribosome"]))
    long = label_width_pt(
        pd.Series(["protein transport", "positive regulation of gene expression"])
    )
    assert long > short > 0


def test_rotated_label_height_pt_collapses_to_one_line_at_zero_degrees():
    """Unrotated, a label band is one tick line tall; the sine/cosine terms must not swap."""
    names = pd.Series(["kmeans_merge9_manual"])
    assert rotated_label_height_pt(names, 0) == pytest.approx(
        rotated_label_height_pt(pd.Series(["x"]), 0)
    )
    assert rotated_label_height_pt(names, 90) == pytest.approx(label_width_pt(names))
