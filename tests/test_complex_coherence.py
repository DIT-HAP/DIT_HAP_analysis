"""Tests for complex coherence algorithm (Weiszfeld geometric median + permutation test)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.src.complex.coherence import (
    geometric_median,
    coherence_metrics,
    compute_distance_zscore,
    EPSILON,
)


def test_epsilon_guard_value():
    """Zero-distance epsilon is exactly 1e-5 (source notebook quirk)."""
    assert EPSILON == 1e-5


def test_geometric_median_single_point():
    """geometric_median of a single point is that point itself."""
    points = np.array([[1.0, 2.0]])
    gm = geometric_median(points)
    np.testing.assert_allclose(gm, [1.0, 2.0], atol=1e-6)


def test_geometric_median_collinear_symmetric():
    """geometric_median of symmetric points converges to centroid."""
    points = np.array([[-1.0, 0.0], [1.0, 0.0], [0.0, 0.0]])
    gm = geometric_median(points)
    np.testing.assert_allclose(gm, [0.0, 0.0], atol=1e-4)


def test_geometric_median_uses_component_wise_median_init():
    """Initialization from component-wise median (not mean)."""
    # With points that have outliers, median init is more robust than mean init.
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 0.0], [100.0, 0.0]])
    gm = geometric_median(points)
    # Should converge near the cluster of 3, not pulled to outlier
    assert gm[0] < 2.0


def test_coherence_metrics_returns_required_keys():
    """coherence_metrics returns centroid_x, centroid_y + 6 distance stats."""
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 0.866], [0.5, -0.866]])
    result = coherence_metrics(points)
    required = {"centroid_x", "centroid_y", "median_distance", "mean_distance",
                "std_distance", "min_distance", "max_distance", "mpd"}
    assert required.issubset(set(result.keys()))


def test_compute_distance_zscore_returns_zscore_and_pvalue():
    """compute_distance_zscore returns observed_mpd, z_score, p_value, n_permutations."""
    rng = np.random.default_rng(42)
    all_points = rng.standard_normal((100, 2))
    complex_indices = list(range(5))  # tight cluster
    result = compute_distance_zscore(all_points, complex_indices, n_permutations=100, random_state=42)
    assert "observed_mpd" in result
    assert "z_score" in result
    assert "p_value" in result
    assert result["n_permutations"] == 100


def test_compute_distance_zscore_tight_cluster_has_low_zscore():
    """A tight cluster should have a lower MPD than random (negative z-score)."""
    rng = np.random.default_rng(42)
    # Background: spread out
    all_points = rng.standard_normal((200, 2)) * 5
    # Tight complex: 8 points near origin
    all_points[:8] = rng.standard_normal((8, 2)) * 0.1
    result = compute_distance_zscore(all_points, list(range(8)), n_permutations=500, random_state=42)
    assert result["z_score"] < 0  # observed MPD < permutation mean
