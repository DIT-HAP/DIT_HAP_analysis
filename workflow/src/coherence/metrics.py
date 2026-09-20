"""Geometric coherence metrics for a set of genes in normalized DR-DL space.

Ported faithfully from DIT_HAP_pipeline complex_analysis.ipynb (cells 3, 30, 34).
A gene set is "coherent" when its members sit tighter in (DR, DL) space than a
random gene set of equal size — quantified as a permutation z-score (negative =
tighter than random = coherent). The primary axis is median_pairwise_distance.

Normalization: DR_NORM_MAX / DL_NORM_MAX are plain divisors, so normalized_DR =
DR and normalized_DL = DL / 10. Every metric here is a Euclidean distance in
that space, and a reflection of one axis leaves all distances — and therefore
every z-score — unchanged. Upstream flipped the DR sign on 2026-09-17 (negative
DR is now the depleted end), which mirrors the space without moving the
coherence results; only the plots' DR axis orientation changes.

This module is the single source of truth for both the z-scores and the raw
descriptive statistics: `compute_distance_zscores` handles every dispersion
method in one permutation pass, and `coherence_metrics` returns the per-group
statistic block the workflow's metrics tables are built from.
"""
# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from __future__ import annotations

# 2. Data Processing Imports
import numpy as np

# 3. Third-party Imports
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist, pdist


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
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

# Which `pairwise_distance` reducer backs each pairwise entry of ZSCORE_METHODS.
# Used to decide whether a draw's distance matrix can be shared across methods.
_PAIRWISE_REDUCERS = {
    "mean_pairwise_distance": "mean",
    "median_pairwise_distance": "median",
    "max_pairwise_distance": "max",
}


# =============================================================================
# CORE LOGIC
# =============================================================================
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


def distance_to_centroid(
    X: np.ndarray, centroid: np.ndarray | None = None, method: str = "median"
) -> float | np.ndarray:
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


def pairwise_distance(
    X: np.ndarray, method: str = "median", k_nn: int = 3, pw: np.ndarray | None = None
) -> float | np.ndarray:
    """Pairwise distances between points, reduced by `method`."""
    # `pw` is an optional precomputed `pdist(X)`. Callers that need several pairwise
    # statistics of the same point set pass it so the O(n^2) distance matrix is
    # built once instead of per statistic.
    if pw is None:
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


def _observed_distance(X: np.ndarray, method: str, pw: np.ndarray | None = None) -> float:
    """Observed dispersion of X under `method`; `pw` = a precomputed pdist(X), if any."""
    if method in _PAIRWISE_REDUCERS:
        return pairwise_distance(X, method=_PAIRWISE_REDUCERS[method], pw=pw)
    if method == "mean_distance_to_centroid":
        return distance_to_centroid(X, method="mean")
    if method == "median_distance_to_centroid":
        return distance_to_centroid(X, method="median")
    if method == "mean_knn_distance":
        return average_knn_distance(X, k=3, method="mean")
    raise ValueError(f"Invalid method {method!r}")


def compute_distance_zscores(
    X: np.ndarray,
    bg: np.ndarray,
    methods: list[str] | tuple[str, ...],
    n_permutations: int = DEFAULT_N_PERMUTATIONS,
    random_state: int | None = DEFAULT_RANDOM_STATE,
) -> dict[str, tuple[float, float]]:
    """Permutation z-scores for SEVERAL dispersion methods in one permutation pass."""
    # Returns {method: (z_score, p_value)}. Negative z = tighter than random =
    # coherent. This is a one-sided test for coherence (dispersion tighter than
    # random). The p-value uses the standard add-one estimator
    # p = (#{null <= observed} + 1) / (n_permutations + 1): the +1 counts the
    # observed statistic itself as one of its own permutations, which (a) keeps p
    # strictly positive — a permutation p can never truly be 0, its resolution
    # floor is 1/(n_permutations + 1) — and (b) makes the estimator unbiased, so it
    # survives downstream multiple-testing (FDR) correction instead of collapsing
    # spuriously exact zeros. A method gets (0.0, 1.0) for n<=1 or a zero-variance
    # null.
    #
    # Every method shares the SAME null draws: the rng is seeded once and each
    # permutation picks one index set, which is then scored by every method. That
    # is exactly what scoring the methods one at a time would produce (same seed,
    # same draw sequence), so results are unaffected — only the repeated draws and
    # the repeated pairwise distance matrix go away. When any requested method is
    # pairwise-reduced, each draw's `pdist` is computed once and shared.
    methods = list(dict.fromkeys(methods))
    n_samples = X.shape[0]
    if n_samples <= 1:
        return {method: (0.0, 1.0) for method in methods}

    # One pdist for the whole draw is only worth building if something consumes it.
    shares_pw = any(method in _PAIRWISE_REDUCERS for method in methods)

    rng = np.random.default_rng(random_state)
    idx = np.arange(bg.shape[0])
    observed = {
        method: _observed_distance(X, method, pdist(X) if shares_pw else None)
        for method in methods
    }
    nulls = {method: np.empty(n_permutations) for method in methods}
    for i in range(n_permutations):
        pick = rng.choice(idx, size=n_samples, replace=False)
        draw = bg[pick]
        pw = pdist(draw) if shares_pw else None
        for method in methods:
            nulls[method][i] = _observed_distance(draw, method, pw)

    results: dict[str, tuple[float, float]] = {}
    for method in methods:
        null = nulls[method]
        p_value = float((np.sum(null <= observed[method]) + 1) / (n_permutations + 1))
        std = null.std()
        if std == 0:
            results[method] = (0.0, 1.0)
        else:
            results[method] = (float((observed[method] - null.mean()) / std), p_value)
    return results


def compute_distance_zscore(
    X: np.ndarray,
    bg: np.ndarray,
    method: str,
    n_permutations: int = DEFAULT_N_PERMUTATIONS,
    random_state: int | None = DEFAULT_RANDOM_STATE,
) -> tuple[float, float]:
    """Permutation z-score of a single dispersion metric. See `compute_distance_zscores`."""
    return compute_distance_zscores(X, bg, [method], n_permutations, random_state)[method]


def coherence_metrics(points: np.ndarray) -> dict:
    """Geometric-median location + descriptive stats of all pairwise L2 distances."""
    # The per-group statistic block that `compute_coherence.py` writes into the
    # metrics table. The keys are the column names of that table, so they name the
    # estimand exactly: `geom_median_*` is the GEOMETRIC MEDIAN (Weiszfeld), not a
    # centroid (a centroid is a mean, a different point), and every `*_pairwise_*`
    # column is a reduction over the pdist matrix. The five `*_pairwise_distance`
    # keys share their spelling with `ZSCORE_METHODS`, so a method's observed value
    # and its z/p columns are all derivable from the method name.
    #
    # A single-point group has no pairwise distances; every distance statistic is
    # 0.0 and only the location is meaningful.
    points = np.asarray(points, dtype=float)
    centroid = geometric_median(points)

    pairwise = pdist(points)
    if pairwise.size == 0:
        median_d = mean_d = std_d = min_d = max_d = 0.0
    else:
        median_d = float(np.median(pairwise))
        mean_d = float(np.mean(pairwise))
        std_d = float(np.std(pairwise))
        min_d = float(np.min(pairwise))
        max_d = float(np.max(pairwise))

    return {
        "geom_median_DR": float(centroid[0]),
        "geom_median_DL": float(centroid[1]),
        "median_pairwise_distance": median_d,
        "mean_pairwise_distance": mean_d,
        "std_pairwise_distance": std_d,
        "min_pairwise_distance": min_d,
        "max_pairwise_distance": max_d,
    }
