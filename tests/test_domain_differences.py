"""Tests for intra-gene DR heterogeneity (domain-difference) candidate stats.

The source notebook genes_with_domain_differences.ipynb is a visualization
notebook: it selects genes with gene-level DR > 0.15 and scatter-plots each
in-gene insertion (Residue_affected vs DR). It computes NO per-gene stats
table. compute_domain_candidate_stats is the deterministic distillation of
that visualization into per-gene positional-spread statistics — the pieces
below are the unit-tested primitives it is built from.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.scripts.domain_differences.compute_domain_stats import (
    DR_THRESHOLD,
    IN_GENE_FILTER,
    _normalize_legacy_metrics,
    compute_insertion_fraction,
    filter_high_dr_genes,
    compute_domain_candidate_stats,
    resolve_duplicate_annotations,
)


# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------
def test_dr_threshold_constant():
    """DR_THRESHOLD == 0.15 is the exact gene-selection cutoff from the notebook."""
    assert DR_THRESHOLD == 0.15


def test_in_gene_filter_matches_repo_convention():
    """The in-gene filter string is byte-faithful to the notebook / coverage.smk quirk."""
    assert IN_GENE_FILTER == "Type != 'Intergenic region' and Distance_to_stop_codon > 4"


# ---------------------------------------------------------------------------
# compute_insertion_fraction
# ---------------------------------------------------------------------------
def test_compute_insertion_fraction_basic():
    """insertion_fraction = Distance_to_start_codon / (start + stop)."""
    df = pd.DataFrame({
        "Distance_to_start_codon": [0, 30, 100],
        "Distance_to_stop_codon": [100, 70, 0],
    })
    out = compute_insertion_fraction(df)
    assert out["insertion_fraction"].tolist() == pytest.approx([0.0, 0.3, 1.0])


def test_compute_insertion_fraction_clamped_to_unit_interval():
    """Out-of-range fractions clamp to [0, 1] (negative -> 0, >1 -> 1)."""
    df = pd.DataFrame({
        "Distance_to_start_codon": [-10, 150],
        "Distance_to_stop_codon": [110, -50],  # denom=100 -> -0.1 and 1.5
    })
    out = compute_insertion_fraction(df)
    assert out["insertion_fraction"].tolist() == pytest.approx([0.0, 1.0])


def test_compute_insertion_fraction_zero_denominator_is_nan():
    """A zero-length span (both distances 0) yields NaN rather than inf."""
    df = pd.DataFrame({
        "Distance_to_start_codon": [0],
        "Distance_to_stop_codon": [0],
    })
    out = compute_insertion_fraction(df)
    assert out["insertion_fraction"].isna().all()


def test_compute_insertion_fraction_does_not_mutate_input():
    """Input frame is left untouched (function returns a copy)."""
    df = pd.DataFrame({
        "Distance_to_start_codon": [10],
        "Distance_to_stop_codon": [90],
    })
    _ = compute_insertion_fraction(df)
    assert "insertion_fraction" not in df.columns


# ---------------------------------------------------------------------------
# filter_high_dr_genes
# ---------------------------------------------------------------------------
def test_filter_high_dr_genes_keeps_strictly_above_threshold():
    """Keeps genes with gene-level DR strictly > 0.15; NaN DR is dropped."""
    genes = pd.DataFrame({
        "Systematic ID": ["g1", "g2", "g3", "g4", "g5"],
        "DR": [0.10, 0.15, 0.16, 0.50, np.nan],
    })
    out = filter_high_dr_genes(genes)
    assert out["Systematic ID"].tolist() == ["g3", "g4"]


def test_filter_high_dr_genes_boundary_excluded():
    """DR exactly equal to the threshold is excluded (strict >)."""
    genes = pd.DataFrame({"Systematic ID": ["g1"], "DR": [DR_THRESHOLD]})
    assert len(filter_high_dr_genes(genes)) == 0


# ---------------------------------------------------------------------------
# compute_domain_candidate_stats
# ---------------------------------------------------------------------------
def _make_inputs():
    in_gene = pd.DataFrame({
        "Systematic ID": (
            ["gA", "gA"] + ["gB", "gB", "gB"] + ["gC"] + ["gD"]
        ),
        "insertion_fraction": [0.1, 0.9, 0.4, 0.5, 0.6, 0.5, 0.2],
    })
    high_dr = pd.DataFrame({
        "Systematic ID": ["gA", "gB", "gC"],
        "Name": ["a1", "b1", "c1"],
        "DR": [0.50, 0.20, 0.90],
    })
    return in_gene, high_dr


def test_compute_domain_candidate_stats_columns_and_rows():
    """Output has the expected schema and one row per high-DR gene present."""
    in_gene, high_dr = _make_inputs()
    stats = compute_domain_candidate_stats(in_gene, high_dr)
    assert list(stats.columns) == [
        "Systematic ID", "Name", "n_insertions",
        "mean_insertion_fraction", "std_insertion_fraction", "gene_DR",
    ]
    # gD is not a high-DR gene -> dropped by the inner join.
    assert set(stats["Systematic ID"]) == {"gA", "gB", "gC"}


def test_compute_domain_candidate_stats_values():
    """n_insertions, mean, std, and gene_DR are computed per gene."""
    in_gene, high_dr = _make_inputs()
    stats = compute_domain_candidate_stats(in_gene, high_dr).set_index("Systematic ID")

    assert stats.loc["gA", "n_insertions"] == 2
    assert stats.loc["gB", "n_insertions"] == 3
    assert stats.loc["gC", "n_insertions"] == 1
    assert stats.loc["gA", "mean_insertion_fraction"] == pytest.approx(0.5)
    assert stats.loc["gB", "mean_insertion_fraction"] == pytest.approx(0.5)
    # sample std (ddof=1): gA=sqrt(0.32)=0.5657, gB=0.1, gC single->NaN
    assert stats.loc["gA", "std_insertion_fraction"] == pytest.approx(np.sqrt(0.32))
    assert stats.loc["gB", "std_insertion_fraction"] == pytest.approx(0.1)
    assert np.isnan(stats.loc["gC", "std_insertion_fraction"])
    assert stats.loc["gA", "gene_DR"] == pytest.approx(0.50)
    assert stats.loc["gC", "Name"] == "c1"


def test_compute_domain_candidate_stats_sorted_by_std_desc():
    """Rows are sorted by std_insertion_fraction descending, NaN last."""
    in_gene, high_dr = _make_inputs()
    stats = compute_domain_candidate_stats(in_gene, high_dr)
    # gA (0.566) > gB (0.1) > gC (NaN, last)
    assert stats["Systematic ID"].tolist() == ["gA", "gB", "gC"]


def test_compute_domain_candidate_stats_dedups_gene_key():
    """A duplicated gene key must NOT cross-product (one row, un-inflated count).

    Guards the drop_duplicates('Systematic ID') hardening: real gene-level data
    has unique keys, but a duplicate would otherwise turn both merges into a
    cross-product (e.g. gA's 2 insertions reported as n_insertions=4, two rows).
    """
    in_gene, high_dr = _make_inputs()
    high_dr_dup = pd.concat([high_dr, high_dr.iloc[[0]]], ignore_index=True)  # duplicate gA
    stats = compute_domain_candidate_stats(in_gene, high_dr_dup).set_index("Systematic ID")
    assert (stats.index == ["gA", "gB", "gC"]).all()  # gA appears exactly once
    assert stats.loc["gA", "n_insertions"] == 2       # not 4


# ---------------------------------------------------------------------------
# resolve_duplicate_annotations
# ---------------------------------------------------------------------------
def _dup_indexed_annotations():
    """Two coordinates, the first duplicated with one in-gene + one intergenic row."""
    idx = pd.MultiIndex.from_tuples(
        [("I", 100, "+", "T"), ("I", 100, "+", "T"), ("I", 200, "-", "T")],
        names=["Chr", "Coordinate", "Strand", "Target"],
    )
    return pd.DataFrame({
        "Type": ["Intergenic region", "CDS", "CDS"],
        "Distance_to_stop_codon": [0, 50, 80],
        "Systematic ID": ["gX", "gX", "gY"],
    }, index=idx)


def test_resolve_duplicate_annotations_keeps_in_gene_pass():
    """Among duplicate-indexed rows, the IN_GENE_FILTER-passing row wins."""
    out = resolve_duplicate_annotations(_dup_indexed_annotations())
    assert not out.index.duplicated().any()  # collapsed to unique index
    # The (I,100,+,T) row kept must be the CDS one (passes filter), not intergenic.
    kept = out.loc[[("I", 100, "+", "T")]].iloc[0]
    assert kept["Type"] == "CDS"


def test_resolve_duplicate_annotations_noop_without_duplicates():
    """No duplicate index -> frame returned unchanged."""
    df = _dup_indexed_annotations().iloc[[1, 2]]  # drop the duplicate -> unique index
    out = resolve_duplicate_annotations(df)
    pd.testing.assert_frame_equal(out, df)


# ---------------------------------------------------------------------------
# legacy um/lam -> DR/DL normalization
# ---------------------------------------------------------------------------
def test_normalize_legacy_metrics_renames_um_lam():
    """Legacy um/lam headers are renamed to DR/DL on load."""
    df = pd.DataFrame({"Systematic ID": ["g1"], "um": [0.5], "lam": [3.0]})
    out = _normalize_legacy_metrics(df)
    assert "DR" in out.columns and "DL" in out.columns
    assert "um" not in out.columns and "lam" not in out.columns
    assert out["DR"].iloc[0] == 0.5


def test_normalize_legacy_metrics_idempotent_when_already_dr_dl():
    """Files already shipping DR/DL are left unchanged (no clobber)."""
    df = pd.DataFrame({"Systematic ID": ["g1"], "DR": [0.5], "DL": [3.0]})
    out = _normalize_legacy_metrics(df)
    assert out["DR"].iloc[0] == 0.5 and out["DL"].iloc[0] == 3.0
