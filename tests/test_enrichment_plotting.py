"""Tests for the pure data preparation behind the enrichment figures."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.src.enrichment.plotting import (
    dotplot_matrix,
    dotplot_terms,
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


def test_dotplot_terms_orders_by_namespace_then_significance():
    """BP rows come before CC rows; within a namespace the strongest term leads."""
    terms = dotplot_terms(_two_cluster_frame(), top_terms_per_cluster=2)
    assert list(terms["term_id"]) == ["GO:2", "GO:1", "GO:4"]


def test_dotplot_matrix_maps_rows_and_colours():
    """Every selected row carries its term row index and a -log10(FDR) value."""
    enrichment = _two_cluster_frame()
    terms = dotplot_terms(enrichment, top_terms_per_cluster=2)
    matrix = dotplot_matrix(enrichment, terms)

    assert len(matrix) == 4  # (cluster, term) pairs kept from the selected terms
    assert matrix["row"].max() == len(terms) - 1
    strongest = matrix.loc[matrix["term_id"] == "GO:2", "neg_log10_fdr"].max()
    assert strongest == pytest.approx(6.0)


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
        {"first": _two_cluster_frame(), "second": _enrichment([[1, "BP", "GO:1", "a", 0.5, 0.1]])},
        fdr_threshold=1e-2,
    )
    assert counts.loc[1, "first"] == 2
    assert counts.loc[2, "first"] == 2
    assert counts.loc[1, "second"] == 0


def test_term_display_names_maps_every_seen_term():
    """The name map covers terms from both variants."""
    names = term_display_names({"first": _two_cluster_frame(), "second": _two_cluster_frame()})
    assert names["GO:4"] == "component D"
