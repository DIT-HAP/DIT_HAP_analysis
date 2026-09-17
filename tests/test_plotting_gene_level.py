"""Tests for workflow/src/plotting/gene_level.py."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd

from workflow.src.plotting.gene_level import (
    CLUSTER_COLORS,
    MULTI_COLORS,
    RAW_VALUE_COLS,
    FITTED_VALUE_COLS,
    visualize_cluster_on_feature_space,
    plot_groups_on_feature_space,
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
