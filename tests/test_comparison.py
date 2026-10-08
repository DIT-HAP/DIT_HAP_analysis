"""Tests for large-scale study comparison core computations."""

import numpy as np
import pandas as pd
import pytest

from comparison.core import (
    CLIP_UPPER,
    DENSITY_COLUMNS,
    DISTRIBUTION_SCALES,
    MIN_PAIRS_FOR_CORRELATION,
    STATS_COLUMNS,
    build_distribution_frame,
    build_fitness_table,
    clip_density_columns,
    cluster_column_order,
    compute_correlation_stats,
    compute_correlations,
    correlation_matrix,
    plot_correlation_heatmap,
    plot_fitness_distributions,
    plot_pairwise_scatter,
    rename_metrics_for_comparison,
    select_fitness_columns,
)


def test_clip_upper_constant():
    """clip(upper=200) is the exact value from source notebook."""
    assert CLIP_UPPER == 200


def test_density_columns_include_required_names():
    """Integration density, ipkm, uipkm columns must be clipped."""
    required = {"Integration density, in-vivo (integrations/kb/million inserts)", "ipkm", "uipkm"}
    assert required.issubset(set(DENSITY_COLUMNS))


def test_clip_density_columns_caps_at_200():
    """Values above 200 are clipped to exactly 200."""
    df = pd.DataFrame({
        "Integration density, in-vivo (integrations/kb/million inserts)": [50.0, 250.0, 200.0],
        "ipkm": [100.0, 300.0, 199.0],
        "uipkm": [10.0, 201.0, 5.0],
        "other_col": [1000.0, 2000.0, 3000.0],
    })
    result = clip_density_columns(df)
    assert result["Integration density, in-vivo (integrations/kb/million inserts)"].max() == 200.0
    assert result["ipkm"].max() == 200.0
    assert result["uipkm"].max() == 200.0
    # other_col untouched
    assert result["other_col"].max() == 3000.0


def test_compute_correlations_returns_both_coefficients():
    """compute_correlations returns (pearson r, pearson p, spearman rho, spearman p)."""
    x = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0])
    y = pd.Series([1.1, 1.9, 3.1, 3.9, 5.1])
    r, p, rho, rho_p = compute_correlations(x, y)
    assert abs(r - 1.0) < 0.05
    assert rho > 0.99
    assert p < 0.05


def test_compute_correlations_ignores_nan_pairs():
    """NaN in either column → drop pair before correlation."""
    x = pd.Series([1.0, 2.0, np.nan, 4.0])
    y = pd.Series([1.0, 2.0, 3.0, np.nan])
    r, p, rho, _ = compute_correlations(x, y)
    # Only (1,1) and (2,2) survive — perfect positive correlation
    assert abs(r - 1.0) < 0.01
    assert rho == pytest.approx(1.0)


def test_rename_metrics_requires_both_source_columns():
    """A missing metric column in the annotation reference raises loudly."""
    ref = pd.DataFrame({"HD_DIT_HAP_DR": [1.0]}, index=["SPX1"])
    with pytest.raises(KeyError, match="gRNA_DR"):
        rename_metrics_for_comparison(ref)


def test_build_fitness_table_joins_metrics_on_gene_id():
    """Features spine left-joins both metrics; unmatched genes get NaN metrics."""
    ref = pd.DataFrame(
        {"HD_DIT_HAP_DR": [-1.0], "gRNA_DR": [0.5]},
        index=pd.Index(["SPAC1"], name="gene_systematic_id"),
    )
    features = pd.DataFrame({
        "gene_systematic_id": ["SPAC1", "SPAC2"],
        "Barseq_from_dulab": [0.1, 0.2],
    })
    table = build_fitness_table(ref, features)
    assert len(table) == 2
    assert table.loc[table["gene_systematic_id"] == "SPAC1", "DIT-HAP DR"].item() == -1.0
    assert table.loc[table["gene_systematic_id"] == "SPAC1", "gRNA DR"].item() == 0.5
    # SPAC2 not in the reference -> NaN metrics, features kept
    assert np.isnan(table.loc[table["gene_systematic_id"] == "SPAC2", "DIT-HAP DR"].item())
    assert table.loc[table["gene_systematic_id"] == "SPAC2", "Barseq_from_dulab"].item() == 0.2


def test_compute_correlation_stats_skips_degenerate_and_sparse_pairs():
    """Constant columns and below-threshold overlaps are skipped, not emitted."""
    df = pd.DataFrame({
        "a": [1.0, 2, 3, 4, 5, 6],
        "b": [1.1, 2.2, 3.3, 4.4, 5.5, 6.6],
        "c": [1.0] * 6,  # constant
        "d": [np.nan] * 6,  # empty
    })
    stats = compute_correlation_stats(df, ["a", "b", "c", "d"])
    assert list(stats["pair"]) == ["a vs b"]
    assert set(stats.columns) == set(STATS_COLUMNS)


def test_compute_correlation_stats_bh_correction_monotone():
    """BH-corrected p-values are >= the raw ones and monotone in rank."""
    x = pd.Series(np.arange(20, dtype=float))
    y1 = x + 0.01  # strong
    y2 = x * -1.0  # strong negative
    y3 = pd.Series(np.random.default_rng(0).normal(size=20))  # weak
    df = pd.DataFrame({"x": x, "y1": y1, "y2": y2, "y3": y3})
    stats = compute_correlation_stats(df, ["x", "y1", "y2", "y3"])
    assert (stats["p_fdr"] >= stats["p_pearson"] - 1e-12).all()
    assert (stats["p_spearman_fdr"] >= stats["p_spearman"] - 1e-12).all()
    # Weak pair gets a larger corrected p than the strong pairs
    weak = stats.loc[stats["pair"] == "x vs y3", "p_fdr"].item()
    strong = stats.loc[stats["pair"] == "x vs y1", "p_fdr"].item()
    assert weak >= strong


def test_select_fitness_columns_defensive():
    """Only columns present with >= MIN_PAIRS non-NaN values survive selection."""
    df = pd.DataFrame({
        "Barseq_from_dulab": [1.0, 2.0, np.nan, 3.0],
        "Barseq_from_koch": [np.nan] * 4,
        "DIT-HAP DR": [1.0, 2.0, 3.0, 4.0],
        "gRNA DR": [1.0, 2.0, 3.0, 4.0],
    })
    columns = select_fitness_columns(df)
    assert "Barseq_from_dulab" in columns
    assert "Barseq_from_koch" not in columns
    assert "DIT-HAP DR" in columns and "gRNA DR" in columns


def test_min_pairs_constant():
    """The pairwise threshold matches the notebook's minimal-overlap guard."""
    assert MIN_PAIRS_FOR_CORRELATION == 3


def _two_column_stats() -> pd.DataFrame:
    """Minimal long-form stats over three columns, as compute_correlation_stats emits."""
    return pd.DataFrame({
        "col_x": ["a", "a", "b"],
        "col_y": ["b", "c", "c"],
        "pair": ["a vs b", "a vs c", "b vs c"],
        "n": [10, 10, 10],
        "r_pearson": [0.5, -0.2, 0.9],
        "p_pearson": [0.1, 0.3, 0.01],
        "p_fdr": [0.1, 0.3, 0.01],
        "rho_spearman": [0.6, -0.3, 0.8],
        "p_spearman": [0.05, 0.2, 0.02],
        "p_spearman_fdr": [0.05, 0.2, 0.02],
    })


def test_correlation_matrix_is_symmetric_with_unit_diagonal():
    """Every coefficient lands in both triangles; the diagonal stays 1."""
    matrix = correlation_matrix(_two_column_stats(), ["a", "b", "c"], "r_pearson")
    assert list(matrix.index) == ["a", "b", "c"]
    assert (np.diag(matrix.to_numpy()) == 1).all()
    assert matrix.loc["a", "b"] == matrix.loc["b", "a"] == 0.5
    assert matrix.loc["b", "c"] == 0.9


def test_correlation_matrix_leaves_dropped_pairs_at_zero():
    """A pair the stats TSV does not carry stays 0 rather than raising or NaN-ing."""
    stats = _two_column_stats().iloc[:1]  # only a vs b survives
    matrix = correlation_matrix(stats, ["a", "b", "c"], "r_pearson")
    assert matrix.loc["a", "c"] == 0.0
    assert (np.diag(matrix.to_numpy()) == 1).all()


def test_cluster_column_order_places_similar_columns_together():
    """The order groups the strongly correlated columns and returns a permutation."""
    order = cluster_column_order(_two_column_stats(), ["a", "b", "c"])
    assert sorted(order) == ["a", "b", "c"]
    # b and c correlate at 0.9, against -0.2 / 0.5 for their pairs with a: they
    # cluster first, so a cannot sit between them.
    assert abs(order.index("b") - order.index("c")) == 1


def test_plot_pairwise_scatter_smoke(tmp_path):
    """plot_pairwise_scatter renders the n x n matrix as a single PDF + review PNG."""
    rng = np.random.default_rng(0)
    df = pd.DataFrame({column: rng.normal(size=60) for column in "abcd"})
    columns = ["a", "b", "c", "d"]
    stem = tmp_path / "pairwise_fitness_comparison"
    plot_pairwise_scatter(df, columns, stem, order=["d", "b", "a", "c"])
    assert (tmp_path / "pairwise_fitness_comparison.pdf").exists()
    assert (tmp_path / "pairwise_fitness_comparison.review.png").exists()
    # One figure, not one page per pair.
    assert not list(tmp_path.glob("pairwise_fitness_comparison_p*.pdf"))


def test_plot_correlation_heatmap_smoke(tmp_path):
    """plot_correlation_heatmap renders one two-panel figure for both coefficients."""
    stats = _two_column_stats()
    columns = ["a", "b", "c"]
    output = tmp_path / "correlation_heatmap.pdf"
    plot_correlation_heatmap(stats, columns, output, order=["b", "c", "a"])
    assert output.exists()
    assert output.stat().st_size > 0


def test_build_distribution_frame_pairs_each_column_with_its_log10():
    """Every column appears under both scales; the log10 half drops non-positive values."""
    fitness_table = pd.DataFrame({"a": [1.0, 10.0, 100.0, -5.0, np.nan], "b": [2.0, 4.0, 8.0, 16.0, 32.0]})
    frame = build_distribution_frame(fitness_table, ["a", "b"])
    assert list(dict.fromkeys(frame["phenotype"])) == ["a", "b"]
    assert list(dict.fromkeys(frame["scale"])) == list(DISTRIBUTION_SCALES)

    raw_a = frame[(frame["phenotype"] == "a") & (frame["scale"] == DISTRIBUTION_SCALES[0])]
    log_a = frame[(frame["phenotype"] == "a") & (frame["scale"] == DISTRIBUTION_SCALES[1])]
    assert raw_a["value"].tolist() == [1.0, 10.0, 100.0, -5.0]  # NaN dropped, sign kept
    assert log_a["value"].tolist() == [0.0, 1.0, 2.0]  # -5 and NaN have no log10


def test_plot_fitness_distributions_smoke(tmp_path):
    """plot_fitness_distributions renders the raw | log10 row-per-column QC figure."""
    rng = np.random.default_rng(0)
    fitness_table = pd.DataFrame({column: np.abs(rng.normal(size=40)) for column in "abc"})
    stem = tmp_path / "fitness_distributions_qc"
    plot_fitness_distributions(fitness_table, stem, order=["c", "a", "b"])
    assert (tmp_path / "fitness_distributions_qc.pdf").exists()
    assert (tmp_path / "fitness_distributions_qc.review.png").exists()
