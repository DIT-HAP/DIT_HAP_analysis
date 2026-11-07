"""Tests for workflow/src/figure_render/cluster_scatter.py."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from workflow.src.figure_render.cluster_scatter import (
    build_cluster_grid_figure,
    build_cluster_scatter_figure,
)


def _synthetic(n: int = 120, n_clusters: int = 4) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    return pd.DataFrame(
        {
            "DR": rng.normal(-0.6, 0.3, n),
            "DL": rng.random(n) * 6,
            "cluster": rng.integers(1, n_clusters + 1, n),
        }
    )


def test_build_cluster_scatter_figure_draws_two_panels_and_one_legend():
    """The per-variant page keeps the density + flat-colour panels and a shared legend."""
    figure = build_cluster_scatter_figure(
        _synthetic(), cluster_col="cluster", one_based=True, legend=True
    )

    assert len(figure.axes) == 2
    assert len(figure.legends) == 1
    plt.close(figure)


def test_build_cluster_grid_figure_hides_unused_cells():
    """A 3-variant grid on 2 columns leaves the fourth cell hidden."""
    frames = {"a": _synthetic(), "b": _synthetic(), "c": _synthetic()}

    figure = build_cluster_grid_figure(frames, col_num=2, title="grid")

    assert len([ax for ax in figure.axes if ax.get_visible()]) == 3
    plt.close(figure)
