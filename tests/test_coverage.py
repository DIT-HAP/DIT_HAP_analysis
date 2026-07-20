"""Tests for gene coverage computation logic."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from workflow.scripts.coverage.compute_coverage_stats import (
    IN_GENE_FILTER,
    compute_insertion_coverage,
    compute_gene_coverage,
    compute_essentiality_coverage,
)


def _make_insertion_annotation(n_in_gene=30, n_intergenic=10):
    """Synthetic insertion annotation table with required columns."""
    rows = []
    for i in range(n_in_gene):
        rows.append({"Type": "Coding exon", "Distance_to_stop_codon": 10})
    for i in range(n_intergenic):
        rows.append({"Type": "Intergenic region", "Distance_to_stop_codon": 0})
    # Edge: in-gene but too close to stop codon
    rows.append({"Type": "Coding exon", "Distance_to_stop_codon": 3})
    idx = pd.MultiIndex.from_tuples(
        [(f"I", i * 100, "+", f"g{i}") for i in range(len(rows))],
        names=["Chr", "Coordinate", "Strand", "Gene"],
    )
    return pd.DataFrame(rows, index=idx)


def test_in_gene_filter_constant():
    """Exact filter string is preserved from source notebook (quirk)."""
    assert IN_GENE_FILTER == "Type != 'Intergenic region' and Distance_to_stop_codon > 4"


def test_compute_insertion_coverage_counts():
    """In-gene count = rows passing filter; intergenic = complement."""
    annotation = _make_insertion_annotation(n_in_gene=30, n_intergenic=10)
    # 30 in-gene with Distance_to_stop_codon=10, 1 edge with Distance=3 (fails), 10 intergenic
    result = compute_insertion_coverage(annotation)
    assert result["total"] == 41
    assert result["in_gene"] == 30  # edge case excluded
    assert result["intergenic"] == 11


def test_compute_gene_coverage_counts():
    """covered = DR not NaN; not_covered = DR is NaN."""
    gene_result = pd.DataFrame({
        "Systematic ID": ["g1", "g2", "g3", "g4"],
        "DR": [0.5, None, 0.8, None],
        "DeletionLibrary_essentiality": ["E", "V", "E", "V"],
    })
    result = compute_gene_coverage(gene_result)
    assert result["total"] == 4
    assert result["covered"] == 2
    assert result["not_covered"] == 2


def test_compute_essentiality_coverage_essential():
    """Essential (E) gene coverage split is correct."""
    # 3 essential genes (g0-g2) with 2 covered (DR not-NaN), 3 non-essential
    # genes (g3-g5) with 2 covered. NOTE: fixed from the original migration
    # plan's fixture, which had DR=[0.5, None, 0.8, None, 0.7, None] — that
    # data has only 1 non-null DR in the V group (g3-g5), making
    # non_essential["covered"] == 2 mathematically unsatisfiable regardless
    # of implementation (3 non-null DR values total across all 6 genes, but
    # the two assertions below sum to 4). Corrected here to g4's DR staying
    # non-null and g3 also covered, matching the assertions' intended shape.
    gene_result = pd.DataFrame({
        "Systematic ID": [f"g{i}" for i in range(6)],
        "DR": [0.5, None, 0.8, 0.6, 0.7, None],
        "DeletionLibrary_essentiality": ["E", "E", "E", "V", "V", "V"],
    })
    result = compute_essentiality_coverage(gene_result)
    assert result["essential"]["total"] == 3
    assert result["essential"]["covered"] == 2
    assert result["non_essential"]["total"] == 3
    assert result["non_essential"]["covered"] == 2
