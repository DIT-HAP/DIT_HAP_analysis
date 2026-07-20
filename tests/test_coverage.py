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
    load_gene_level,
    resolve_duplicate_annotations,
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


def test_compute_essentiality_coverage_excludes_not_determined():
    """Genes with essentiality 'Not_determined' land in neither essential nor non_essential.

    Real releases (e.g. HD_DIT_HAP) carry a third essentiality value besides
    'E'/'V'. Splitting on `== 'E'` vs `== 'V'` (matching the source notebook
    and _HIST_ROW_QUERIES) means such genes are excluded from both buckets,
    so essential.total + non_essential.total < len(gene_result).
    """
    gene_result = pd.DataFrame({
        "Systematic ID": ["g1", "g2", "g3", "g4", "g5"],
        "DR": [0.5, 0.6, 0.7, None, 0.9],
        "DeletionLibrary_essentiality": ["E", "V", "Not_determined", "Not_determined", "E"],
    })
    result = compute_essentiality_coverage(gene_result)
    assert result["essential"]["total"] == 2
    assert result["non_essential"]["total"] == 1
    # 2 (essential) + 1 (non_essential) + 2 (Not_determined) == 5 total genes
    assert result["essential"]["total"] + result["non_essential"]["total"] < len(gene_result)


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


def test_load_gene_level_renames_legacy_um_lam(tmp_path):
    """Legacy um/lam headers are renamed to DR/DL."""
    legacy_tsv = tmp_path / "fitting_results.tsv"
    pd.DataFrame({
        "Systematic ID": ["g1", "g2"],
        "um": [0.5, 0.6],
        "lam": [1.0, 2.0],
        "DeletionLibrary_essentiality": ["E", "V"],
    }).to_csv(legacy_tsv, sep="\t", index=False)

    result = load_gene_level(legacy_tsv)
    assert "DR" in result.columns
    assert "DL" in result.columns
    assert "um" not in result.columns
    assert "lam" not in result.columns
    assert list(result["DR"]) == [0.5, 0.6]
    assert list(result["DL"]) == [1.0, 2.0]


def test_load_gene_level_is_idempotent_when_dr_dl_already_present(tmp_path):
    """Rename only triggers when DR/DL aren't already present — a no-op on current-schema files."""
    current_tsv = tmp_path / "fitting_results.tsv"
    pd.DataFrame({
        "Systematic ID": ["g1", "g2"],
        "DR": [0.5, 0.6],
        "DL": [1.0, 2.0],
        "DeletionLibrary_essentiality": ["E", "V"],
    }).to_csv(current_tsv, sep="\t", index=False)

    result = load_gene_level(current_tsv)
    assert list(result.columns) == ["Systematic ID", "DR", "DL", "DeletionLibrary_essentiality"]
    assert list(result["DR"]) == [0.5, 0.6]
    assert list(result["DL"]) == [1.0, 2.0]


def test_resolve_duplicate_annotations_keeps_passing_duplicate():
    """When one duplicate passes IN_GENE_FILTER and one doesn't, the passing one is kept."""
    idx = pd.MultiIndex.from_tuples(
        [("I", 100, "+", "TTAA"), ("I", 100, "+", "TTAA")],
        names=["Chr", "Coordinate", "Strand", "Target"],
    )
    annotations = pd.DataFrame(
        {
            "Type": ["Intergenic region", "Coding gene"],
            "Distance_to_stop_codon": [0, 10],
        },
        index=idx,
    )
    result = resolve_duplicate_annotations(annotations)
    assert len(result) == 1
    # The row that passes IN_GENE_FILTER (Coding gene, Distance=10) is kept.
    assert result.iloc[0]["Type"] == "Coding gene"
    assert result.iloc[0]["Distance_to_stop_codon"] == 10


def test_resolve_duplicate_annotations_deterministic_when_neither_passes():
    """When neither duplicate passes IN_GENE_FILTER, the first row (original file order) is kept deterministically."""
    idx = pd.MultiIndex.from_tuples(
        [("I", 200, "+", "TTAA"), ("I", 200, "+", "TTAA")],
        names=["Chr", "Coordinate", "Strand", "Target"],
    )
    annotations = pd.DataFrame(
        {
            "Type": ["Intergenic region", "Intergenic region"],
            "Distance_to_stop_codon": [0, 0],
            "_marker": ["first", "second"],
        },
        index=idx,
    )
    result = resolve_duplicate_annotations(annotations)
    assert len(result) == 1
    # Neither passes IN_GENE_FILTER, so the stable sort preserves original
    # order and the first row ("first") is kept — deterministic, not
    # dependent on pandas' default (unstable) sort algorithm.
    assert result.iloc[0]["_marker"] == "first"


def test_resolve_duplicate_annotations_no_duplicates_is_noop():
    """No duplicate index values -> annotations pass through unchanged."""
    idx = pd.MultiIndex.from_tuples(
        [("I", 100, "+", "TTAA"), ("I", 200, "+", "TTAA")],
        names=["Chr", "Coordinate", "Strand", "Target"],
    )
    annotations = pd.DataFrame(
        {"Type": ["Coding gene", "Intergenic region"], "Distance_to_stop_codon": [10, 0]},
        index=idx,
    )
    result = resolve_duplicate_annotations(annotations)
    pd.testing.assert_frame_equal(result, annotations)
