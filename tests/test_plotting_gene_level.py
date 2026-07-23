"""Tests for workflow/src/plotting/gene_level.py."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from workflow.src.plotting.gene_level import (
    CLUSTER_COLORS,
    MULTI_COLORS,
    RAW_VALUE_COLS,
    FITTED_VALUE_COLS,
    GRNA_VALUE_COLS,
    visualize_cluster_on_feature_space,
    plot_groups_on_feature_space,
    sigmoid_gompertz,
    plot_gene_depletion_curve,
)


def _synthetic(n: int, n_clusters: int, cluster_col: str) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    return pd.DataFrame(
        {
            "DR": rng.random(n),
            "DL": rng.random(n) * 10,
            cluster_col: rng.integers(1, n_clusters + 1, n),
            "Name": [f"g{i}" for i in range(n)],
        }
    )


def test_value_column_constants():
    """Raw and fitted LFC columns are the 5 DIT-HAP time points."""
    assert RAW_VALUE_COLS == ["YES0", "YES1", "YES2", "YES3", "YES4"]
    assert FITTED_VALUE_COLS == ["YES0_fitted", "YES1_fitted", "YES2_fitted", "YES3_fitted", "YES4_fitted"]


def test_visualize_merged_clusters_two_panels():
    """The merged 1..9 case yields a two-panel figure without palette overrun."""
    df = _synthetic(90, 9, "cluster")
    fig = visualize_cluster_on_feature_space(df, "cluster", show_box=True, cluster_minus_one=True)
    assert len(fig.axes) == 2


def test_visualize_candidate_64_clusters_does_not_overrun_palette():
    """The 64-cluster candidate review must not IndexError on the 10-color palette."""
    df = _synthetic(200, 64, "cluster")
    # cluster_minus_one=False -> idx can reach 63; modulo indexing must keep it safe.
    fig = visualize_cluster_on_feature_space(df, "cluster")
    assert len(fig.axes) == 2


def test_palette_sizes():
    """MULTI_COLORS covers the 64 candidate clusters; CLUSTER_COLORS the 10 merged."""
    assert len(MULTI_COLORS) == 64
    assert len(CLUSTER_COLORS) == 10


def test_plot_groups_on_feature_space_grid_shape():
    """A grid of feature-space subplots has one axis per group (plus filled columns)."""
    df = _synthetic(60, 6, "cluster")
    fig = plot_groups_on_feature_space(df, "cluster", "Name", col_num=3)
    assert len(fig.axes) >= 6


# =============================================================================
# SINGLE-GENE DEPLETION CURVE
# =============================================================================
def test_sigmoid_gompertz_zero_amplitude():
    """A==0 yields a flat zero curve of matching shape."""
    x = np.linspace(0, 13, 50)
    y = sigmoid_gompertz(x, A=0.0, DR=0.5, DL=2.0)
    assert y.shape == x.shape
    assert np.all(y == 0.0)


def test_sigmoid_gompertz_monotone_and_bounded():
    """A non-degenerate curve rises monotonically from ~0 toward the plateau A."""
    x = np.linspace(0, 13, 100)
    A = 6.0
    y = sigmoid_gompertz(x, A=A, DR=0.5, DL=2.0)
    assert np.all(np.diff(y) >= -1e-9)          # non-decreasing
    assert y[0] < 0.5 and y[-1] <= A + 1e-9     # starts low, never exceeds A


def test_sigmoid_gompertz_no_overflow_extreme_dl():
    """A huge DL drives the exponent past the clip bound without overflowing."""
    x = np.linspace(0, 13, 100)
    y = sigmoid_gompertz(x, A=5.0, DR=0.5, DL=1e6)
    assert np.all(np.isfinite(y))


def test_plot_gene_depletion_curve_dit_only():
    """With grna_row=None the panel has 3 lines (fit, slope, DIT-HAP) and no gRNA line."""
    dit_row = pd.Series(
        {"A": 6.0, "DR": 0.5, "DL": 2.0, "YES0": 0.0, "YES1": 1.0, "YES2": 3.0, "YES3": 5.0, "YES4": 6.0}
    )
    fig, ax = plt.subplots()
    plot_gene_depletion_curve(ax, dit_row, grna_row=None, title="myGene")
    labels = [line.get_label() for line in ax.get_lines()]
    assert ax.get_title() == "myGene"
    assert not any("gRNA" == lbl for lbl in labels)
    assert sum(lbl.startswith("DIT_HAP") for lbl in labels) == 1
    plt.close(fig)


def test_plot_gene_depletion_curve_with_grna_overlay():
    """A gRNA row adds a 4th line labeled 'gRNA'."""
    dit_row = pd.Series(
        {"A": 6.0, "DR": 0.5, "DL": 2.0, "YES0": 0.0, "YES1": 1.0, "YES2": 3.0, "YES3": 5.0, "YES4": 6.0}
    )
    grna_row = pd.Series({c: float(i) for i, c in enumerate(GRNA_VALUE_COLS)})
    fig, ax = plt.subplots()
    plot_gene_depletion_curve(ax, dit_row, grna_row=grna_row, title="myGene")
    labels = [line.get_label() for line in ax.get_lines()]
    assert "gRNA" in labels
    plt.close(fig)
