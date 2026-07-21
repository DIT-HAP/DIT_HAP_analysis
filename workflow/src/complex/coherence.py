"""
Complex Coherence Algorithm
===========================

Pure-logic library for the macromolecular-complex coherence analysis: the
Weiszfeld geometric median, pairwise-distance coherence metrics, and a seeded
permutation test for the median pairwise distance (MPD) of a complex's members
against a genome-wide background. Ported from
DIT_HAP_pipeline/workflow/notebooks/complex_analysis.ipynb (section 5.1).

The "points" are genes placed in the 2D DIT-HAP fitness space (min-max
normalized DR/DL); a tight cluster of complex members is more "coherent"
(smaller MPD) than a random draw of the same number of background genes.

This module contains no file I/O and no CLI — just numpy/scipy functions, so
it is cleanly unit-testable (see tests/test_complex_coherence.py).

Usage
-----
    from workflow.src.complex.coherence import (
        geometric_median, coherence_metrics, compute_distance_zscore, EPSILON,
    )

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-20
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
from scipy.spatial.distance import cdist, pdist

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Zero-distance guard for Weiszfeld's inverse-distance weights: any point whose
# distance to the current centre is < EPSILON is clamped to EPSILON so the
# weight 1/distance never blows up. Byte-faithful to the source notebook's
# geometric_median default (epsilon=1e-5).
EPSILON: float = 1e-5

# Convergence tolerance and iteration cap for the Weiszfeld iteration. The
# notebook reused epsilon (1e-5) as the convergence bound; we tighten it to
# 1e-7 for a more precise centre (the migration plan's explicit choice) and add
# an explicit iteration cap so a pathological input can never spin forever.
_CONVERGENCE_TOL: float = 1e-7
_MAX_ITER: int = 1000


# =============================================================================
# GEOMETRIC MEDIAN (Weiszfeld)
# =============================================================================
def geometric_median(points: np.ndarray, epsilon: float = EPSILON) -> np.ndarray:
    """Compute the geometric median of a set of points via Weiszfeld's algorithm.

    The geometric median minimises the sum of Euclidean distances to all
    points and is far more robust to outliers than the arithmetic mean. The
    iteration is initialised from the *component-wise median* (not the mean),
    which keeps the starting point inside the main cluster even when a single
    coordinate is dragged out by an outlier.

    Parameters
    ----------
    points : np.ndarray
        (n_points, n_dims) array of coordinates.
    epsilon : float
        Zero-distance guard (see module-level EPSILON): distances below this
        are clamped so the inverse-distance weights stay finite.

    Returns
    -------
    np.ndarray
        (n_dims,) coordinates of the geometric median.
    """
    points = np.asarray(points, dtype=float)
    if points.shape[0] <= 1:
        # A single point (or empty set with one row) is its own median.
        return points[0].copy()

    # Component-wise median init (robust to per-axis outliers).
    y = np.median(points, axis=0)

    for _ in range(_MAX_ITER):
        # Distance from every point to the current centre; guard zeros so the
        # inverse-distance weights never divide by zero.
        distances = cdist(points, [y]).flatten()
        distances = np.clip(distances, a_min=epsilon, a_max=None)
        weights = 1.0 / distances
        y_next = np.average(points, axis=0, weights=weights)
        if np.linalg.norm(y - y_next) < _CONVERGENCE_TOL:
            y = y_next
            break
        y = y_next

    return y


# =============================================================================
# COHERENCE METRICS
# =============================================================================
def coherence_metrics(points: np.ndarray) -> dict:
    """Summarise the spatial coherence of a set of points.

    Returns the geometric-median centroid coordinates plus summary statistics
    of the full set of pairwise Euclidean (L2) distances. `mpd` (median
    pairwise distance) is the coherence statistic used by the permutation test
    and equals `median_distance`; it is exposed under both names because the
    permutation null is expressed in terms of MPD.

    Degenerate inputs (0 or 1 point) have no pairwise distances; their distance
    statistics are reported as 0.0 (a lone point is maximally coherent), which
    also avoids empty-array warnings from np.median/np.min.

    Parameters
    ----------
    points : np.ndarray
        (n_points, n_dims) array of coordinates.

    Returns
    -------
    dict
        Keys: centroid_x, centroid_y, median_distance, mean_distance,
        std_distance, min_distance, max_distance, mpd.
    """
    points = np.asarray(points, dtype=float)
    centroid = geometric_median(points)

    pairwise = pdist(points)  # condensed vector of all C(n,2) L2 distances
    if pairwise.size == 0:
        # 0 or 1 point: no pairwise distances to summarise.
        median_d = mean_d = std_d = min_d = max_d = 0.0
    else:
        median_d = float(np.median(pairwise))
        mean_d = float(np.mean(pairwise))
        std_d = float(np.std(pairwise))
        min_d = float(np.min(pairwise))
        max_d = float(np.max(pairwise))

    return {
        "centroid_x": float(centroid[0]),
        "centroid_y": float(centroid[1]),
        "median_distance": median_d,
        "mean_distance": mean_d,
        "std_distance": std_d,
        "min_distance": min_d,
        "max_distance": max_d,
        "mpd": median_d,  # median pairwise distance == median_distance
    }


# =============================================================================
# PERMUTATION TEST (median pairwise distance z-score)
# =============================================================================
def _median_pairwise_distance(points: np.ndarray) -> float:
    """Median of all pairwise L2 distances (0.0 when there are none)."""
    pairwise = pdist(np.asarray(points, dtype=float))
    if pairwise.size == 0:
        return 0.0
    return float(np.median(pairwise))


def compute_distance_zscore(
    all_points: np.ndarray,
    complex_indices: list[int],
    n_permutations: int = 1000,
    random_state: int = 42,
) -> dict:
    """Permutation z-score of a complex's median pairwise distance (MPD).

    Compares the observed MPD of the complex's member points against a null
    distribution built by repeatedly drawing `len(complex_indices)` random
    genes (rows) from the genome-wide background (`all_points`) and recomputing
    their MPD. A tight, coherent complex has a smaller MPD than random draws,
    giving a negative z-score and a small (one-sided) p-value.

    The generator is seeded via `np.random.default_rng(random_state)`, so the
    result is deterministic for a fixed seed.

    Parameters
    ----------
    all_points : np.ndarray
        (n_background, n_dims) background point cloud to sample the null from.
    complex_indices : list[int]
        Row indices into `all_points` identifying the complex's member points.
    n_permutations : int
        Number of random draws forming the null distribution.
    random_state : int
        Seed for the permutation RNG.

    Returns
    -------
    dict
        Keys: observed_mpd, z_score, p_value, n_permutations.
    """
    all_points = np.asarray(all_points, dtype=float)
    member_points = all_points[complex_indices]
    n_members = member_points.shape[0]

    observed_mpd = _median_pairwise_distance(member_points)

    # Degenerate: 0 or 1 member -> no meaningful MPD / null distribution.
    if n_members <= 1:
        return {
            "observed_mpd": observed_mpd,
            "z_score": 0.0,
            "p_value": 1.0,
            "n_permutations": n_permutations,
        }

    rng = np.random.default_rng(random_state)
    null_mpds = np.empty(n_permutations, dtype=float)
    for i in range(n_permutations):
        # Draw n_members background rows without replacement (rng.choice on a
        # 2D array samples along the first axis).
        sample = rng.choice(all_points, size=n_members, replace=False)
        null_mpds[i] = _median_pairwise_distance(sample)

    null_mean = float(np.mean(null_mpds))
    null_std = float(np.std(null_mpds))
    # One-sided: fraction of null MPDs at or below the observed (small for a
    # tight cluster).
    p_value = float(np.mean(null_mpds <= observed_mpd))

    # Degenerate null (all permutation MPDs identical): z-score is undefined,
    # report 0.0. Use a tolerance rather than == 0 because np.std returns a
    # tiny non-zero value (~1e-16) for an all-identical array due to float
    # roundoff, which would otherwise produce a spurious huge z-score.
    if null_std < 1e-12:
        z_score = 0.0
    else:
        z_score = (observed_mpd - null_mean) / null_std

    return {
        "observed_mpd": observed_mpd,
        "z_score": z_score,
        "p_value": p_value,
        "n_permutations": n_permutations,
    }
