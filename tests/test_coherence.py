"""Tests for the shared coherence algorithm as used by complex coherence.

Task 6's compute_complex_coherence.py sources the Weiszfeld geometric median and
the seeded median-pairwise-distance permutation test from the canonical
workflow/src/coherence/metrics.py (shared with the themes-A/D verification_complex
scripts). These tests pin the behaviour Task 6 relies on against that module's
API: geometric_median init/convergence, the median_pairwise_distance z-score
sign, and its degenerate-input contract.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pytest

from workflow.src.coherence.metrics import (
    geometric_median,
    coherence_metrics,
    compute_distance_zscore,
    pairwise_distance,
)


def test_geometric_median_single_point():
    """geometric_median of a single point is that point itself."""
    points = np.array([[1.0, 2.0]])
    gm = geometric_median(points)
    np.testing.assert_allclose(gm, [1.0, 2.0], atol=1e-6)


def test_geometric_median_collinear_symmetric():
    """geometric_median of symmetric points converges to the centre."""
    points = np.array([[-1.0, 0.0], [1.0, 0.0], [0.0, 0.0]])
    gm = geometric_median(points)
    np.testing.assert_allclose(gm, [0.0, 0.0], atol=1e-4)


def test_geometric_median_robust_to_outlier():
    """Median-based centre is not dragged out by a single far outlier."""
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 0.0], [100.0, 0.0]])
    gm = geometric_median(points)
    assert gm[0] < 2.0  # stays near the cluster of 3, not pulled toward 100


def test_coherence_metrics_returns_geomedian_and_pairwise_distance_keys():
    """coherence_metrics exposes the geometric median + the pairwise-distance stat block.

    These are the columns compute_coherence.py writes per group, and their names
    are the estimands: `geom_median_*` is the Weiszfeld geometric median (not a
    centroid), and each `*_pairwise_distance` is a reduction over the pdist
    matrix.
    """
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 0.866], [0.5, -0.866]])
    result = coherence_metrics(points)
    for key in ("geom_median_DR", "geom_median_DL", "median_pairwise_distance",
                "mean_pairwise_distance", "std_pairwise_distance",
                "min_pairwise_distance", "max_pairwise_distance"):
        assert key in result
    assert (result["min_pairwise_distance"] <= result["median_pairwise_distance"]
            <= result["max_pairwise_distance"])


def test_compute_distance_zscores_matches_per_method_calls():
    """Merging several methods into ONE permutation pass must not move any number.

    The merged pass seeds the rng once and scores every method on the same draws,
    which is exactly the sequence of draws each single-method call would make on
    its own. Equality is asserted exactly, not approximately — a shared pairwise
    distance matrix must not perturb the arithmetic.
    """
    from workflow.src.coherence.metrics import compute_distance_zscores

    rng = np.random.default_rng(7)
    bg = rng.standard_normal((300, 2))
    for methods in (
        ["median_pairwise_distance", "mean_pairwise_distance"],
        ["mean_knn_distance", "median_pairwise_distance"],
        ["mean_distance_to_centroid", "median_distance_to_centroid"],
    ):
        X = rng.standard_normal((12, 2))
        merged = compute_distance_zscores(X, bg, methods, n_permutations=200, random_state=42)
        for method in methods:
            separate = compute_distance_zscore(
                X, bg, method=method, n_permutations=200, random_state=42
            )
            assert merged[method] == separate, f"{method} moved when merged"


def test_compute_distance_zscores_deduplicates_repeated_methods():
    """Asking for the same method twice scores it once and reports it once."""
    from workflow.src.coherence.metrics import compute_distance_zscores

    rng = np.random.default_rng(11)
    bg = rng.standard_normal((100, 2))
    X = bg[:5]
    merged = compute_distance_zscores(
        X, bg, ["median_pairwise_distance", "median_pairwise_distance"],
        n_permutations=50, random_state=42,
    )
    assert list(merged) == ["median_pairwise_distance"]


def test_compute_distance_zscore_returns_z_and_p_tuple():
    """compute_distance_zscore returns a (z_score, p_value) tuple."""
    rng = np.random.default_rng(42)
    bg = rng.standard_normal((100, 2))
    members = bg[:5]  # a subset of the background
    result = compute_distance_zscore(
        members, bg, method="median_pairwise_distance", n_permutations=100, random_state=42
    )
    assert isinstance(result, tuple) and len(result) == 2
    z, p = result
    assert np.isfinite(z)
    assert 0.0 <= p <= 1.0


def test_compute_distance_zscore_tight_cluster_has_negative_z():
    """A tight complex has a lower MPD than random draws -> negative z-score."""
    rng = np.random.default_rng(42)
    bg = rng.standard_normal((200, 2)) * 5
    bg[:8] = rng.standard_normal((8, 2)) * 0.1  # tight complex near origin
    members = bg[:8]
    z, _ = compute_distance_zscore(
        members, bg, method="median_pairwise_distance", n_permutations=500, random_state=42
    )
    assert z < 0


def test_compute_distance_zscore_single_member_is_degenerate():
    """A 1-member set has no meaningful dispersion -> (0.0, 1.0) per main's contract."""
    rng = np.random.default_rng(0)
    bg = rng.standard_normal((50, 2))
    z, p = compute_distance_zscore(
        bg[:1], bg, method="median_pairwise_distance", n_permutations=100, random_state=42
    )
    assert z == 0.0
    assert p == 1.0


def test_compute_distance_zscore_p_value_never_zero():
    """The add-one estimator keeps p strictly positive at its 1/(n+1) floor.

    A maximally-tight cluster beats every one of the null draws, so #{null<=obs}
    is 0; the add-one term then yields exactly 1/(n_permutations+1), never 0.0
    (a spurious exact zero would poison downstream FDR correction).
    """
    rng = np.random.default_rng(42)
    bg = rng.standard_normal((300, 2)) * 5
    bg[:6] = rng.standard_normal((6, 2)) * 1e-4  # essentially coincident members
    members = bg[:6]
    n_perm = 200
    _, p = compute_distance_zscore(
        members, bg, method="median_pairwise_distance", n_permutations=n_perm, random_state=1
    )
    assert p > 0.0
    assert p == pytest.approx(1.0 / (n_perm + 1))


def test_observed_mpd_equals_median_pairwise_distance():
    """Task 6's observed_mpd (local median pairwise distance) matches the module's."""
    rng = np.random.default_rng(7)
    members = rng.standard_normal((8, 2))
    from scipy.spatial.distance import pdist

    local_mpd = float(np.median(pdist(members)))
    module_mpd = pairwise_distance(members, method="median")
    assert np.isclose(local_mpd, module_mpd)
