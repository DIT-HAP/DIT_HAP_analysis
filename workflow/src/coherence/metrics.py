"""Geometric coherence metrics for a set of genes in normalized DR-DL space.

Ported faithfully from DIT_HAP_pipeline complex_analysis.ipynb (cells 3, 30, 34).
A gene set is "coherent" when its members sit tighter in (DR, DL) space than a
random gene set of equal size — quantified as a permutation z-score (negative =
tighter than random = coherent). The primary axis is median_pairwise_distance.

Normalization (byte-faithful to the notebook): normalized_DR = DR (bounds 0..1),
normalized_DL = DL / 10 (bounds 0..10) — the same DL/10 scaling clustering uses.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial.distance import cdist, pdist
from scipy.spatial import cKDTree

DR_NORM_MAX = 1.0
DL_NORM_MAX = 10.0
DEFAULT_N_PERMUTATIONS = 1000
DEFAULT_RANDOM_STATE = 42
ZSCORE_METHODS = (
    "mean_distance_to_centroid",
    "median_distance_to_centroid",
    "mean_pairwise_distance",
    "median_pairwise_distance",
    "max_pairwise_distance",
    "mean_knn_distance",
)


def normalize_dr_dl(dr: np.ndarray, dl: np.ndarray) -> np.ndarray:
    """Stack DR, DL/10 into an (n, 2) matrix (fixed-bound min-max, min=0)."""
    return np.column_stack([np.asarray(dr, float) / DR_NORM_MAX, np.asarray(dl, float) / DL_NORM_MAX])


def geometric_median(X: np.ndarray, epsilon: float = 1e-5) -> np.ndarray:
    """Weiszfeld's algorithm for the geometric median of points X."""
    y = np.median(X, axis=0)
    while True:
        distances = cdist(X, [y]).flatten()
        distances = np.clip(distances, a_min=epsilon, a_max=None)
        weights = 1.0 / distances
        y_next = np.average(X, axis=0, weights=weights)
        if np.linalg.norm(y - y_next) < epsilon:
            break
        y = y_next
    return y


def average_knn_distance(X: np.ndarray, k: int = 2, method: str = "mean") -> float:
    """Average distance to the k nearest neighbours for each point (KD-tree)."""
    n_samples = X.shape[0]
    if n_samples <= 1:
        return 0.0
    actual_k = min(k, n_samples - 1)
    tree = cKDTree(X)
    distances, _ = tree.query(X, k=actual_k + 1)
    knn = distances[:, 1] if actual_k == 1 else distances[:, 1:]
    return float(np.median(knn)) if method == "median" else float(np.mean(knn))


def distance_to_centroid(X: np.ndarray, centroid: np.ndarray | None = None, method: str = "median"):
    """Distance from each point to the geometric median, reduced by `method`."""
    if centroid is None:
        centroid = geometric_median(X)
    distances = cdist(X, [centroid]).flatten()
    reducers = {
        "mean": lambda: float(np.mean(distances)),
        "median": lambda: float(np.median(distances)),
        "max": lambda: float(np.max(distances)),
        "both": lambda: np.array([np.mean(distances), np.median(distances)]),
        "all": lambda: np.array([np.mean(distances), np.median(distances), np.max(distances)]),
    }
    if method not in reducers:
        raise ValueError(f"Invalid method {method!r}")
    return reducers[method]()


def pairwise_distance(X: np.ndarray, method: str = "median", k_nn: int = 3):
    """Pairwise distances between points, reduced by `method`."""
    pw = pdist(X)
    reducers = {
        "mean": lambda: float(np.mean(pw)),
        "median": lambda: float(np.median(pw)),
        "max": lambda: float(np.max(pw)),
        "knn": lambda: average_knn_distance(X, k=k_nn, method="mean"),
        "both": lambda: np.array([np.mean(pw), np.median(pw)]),
        "all": lambda: np.array([np.mean(pw), np.median(pw), np.max(pw), average_knn_distance(X, k=k_nn)]),
    }
    if method not in reducers:
        raise ValueError(f"Invalid method {method!r}")
    return reducers[method]()


def _observed_distance(X: np.ndarray, method: str) -> float:
    if method == "mean_distance_to_centroid":
        return distance_to_centroid(X, method="mean")
    if method == "median_distance_to_centroid":
        return distance_to_centroid(X, method="median")
    if method == "mean_pairwise_distance":
        return pairwise_distance(X, method="mean")
    if method == "median_pairwise_distance":
        return pairwise_distance(X, method="median")
    if method == "max_pairwise_distance":
        return pairwise_distance(X, method="max")
    if method == "mean_knn_distance":
        return pairwise_distance(X, method="knn")
    raise ValueError(f"Invalid method {method!r}")


def compute_distance_zscore(
    X: np.ndarray,
    bg: np.ndarray,
    method: str,
    n_permutations: int = DEFAULT_N_PERMUTATIONS,
    random_state: int | None = DEFAULT_RANDOM_STATE,
) -> tuple[float, float]:
    """Permutation z-score of a dispersion metric vs random equal-size gene sets.

    Negative z = tighter than random = coherent. p = fraction of permutations
    with dispersion <= observed. Returns (0.0, 1) for n<=1 or zero-variance null.
    """
    rng = np.random.default_rng(random_state)
    n_samples = X.shape[0]
    if n_samples <= 1:
        return 0.0, 1
    observed = _observed_distance(X, method)
    permuted = np.empty(n_permutations)
    idx = np.arange(bg.shape[0])
    for i in range(n_permutations):
        pick = rng.choice(idx, size=n_samples, replace=False)
        permuted[i] = _observed_distance(bg[pick], method)
    std = permuted.std()
    p_value = float(np.mean(permuted <= observed))
    if std == 0:
        return 0.0, 1
    return float((observed - permuted.mean()) / std), p_value


def coherence_metrics(X: np.ndarray) -> dict:
    """All raw dispersion statistics for a gene set (no permutation)."""
    centroid = geometric_median(X)
    mean_dc, median_dc = distance_to_centroid(X, method="both").tolist()
    mean_pw, median_pw, max_pw, mean_knn = pairwise_distance(X, method="all").tolist()
    return {
        "centroid_x": float(centroid[0]),
        "centroid_y": float(centroid[1]),
        "mean_distance_to_centroid": mean_dc,
        "median_distance_to_centroid": median_dc,
        "mean_pairwise_distance": mean_pw,
        "median_pairwise_distance": median_pw,
        "max_pairwise_distance": max_pw,
        "mean_knn_distance": mean_knn,
    }
