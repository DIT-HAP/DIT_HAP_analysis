"""
Clustering grid search (methods x k)
====================================

The DR-DL 2-D clustering comparison behind the redesigned clustering stage:
one method registry, one k range, every partition recorded. The five methods
(KMeans, k-medoids PAM, and the ward/complete/average agglomerative linkages)
all consume the same scaled matrix and differ only in the estimator, so the
module is one dispatch table rather than five near-identical functions.

The k-medoids implementation is a hand-written PAM (build + exact swap) on a
precomputed Euclidean distance matrix. No scikit-learn-extra on this machine,
and PAM is ~45 lines — cheaper than adding a dependency.

Input
-----
- The scaled (DR, DL) matrix pickle from prepare_clustering_data (label
  computation happens on the matrix; labelling of genes happens on its index).

Output
------
- Helper functions returning label arrays / metric dicts; drivers persist them.

Usage
-----
    from clustering.grid_search import fit_method, METHODS
    labels = fit_method("kmedoids", scaled.values, k=9, random_state=42)

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-10-08
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from collections.abc import Callable

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from scipy.cluster.hierarchy import fcluster, linkage as scipy_linkage
from sklearn.cluster import KMeans
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_samples,
    silhouette_score,
)

# 4. Local Imports
from clustering.candidates import renumber_by_dr

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# The five methods, in the fixed order the QC tables and figure report them.
# HC_<linkage> names the family: they all run the same agglomerative skeleton
# and differ only in the cluster-to-cluster distance. complete/average keep the
# old pipeline's pairing with whichever metric it used (complete=cityblock).
METHODS: tuple[str, ...] = ("kmeans", "kmedoids", "HC_ward", "HC_complete", "HC_average")

HC_METRICS: dict[str, dict[str, str]] = {
    "HC_ward": {"linkage": "ward", "metric": "euclidean"},
    "HC_complete": {"linkage": "complete", "metric": "cityblock"},
    "HC_average": {"linkage": "average", "metric": "euclidean"},
}


# =============================================================================
# CORE LOGIC — estimators
# =============================================================================
def cluster_kmeans(coordinates: np.ndarray, k: int, random_state: int) -> np.ndarray:
    """KMeans partition; n_init=10 because a single draw can land two centroids in the DL=0 spike."""
    return KMeans(n_clusters=k, n_init=10, random_state=random_state).fit_predict(coordinates)


def _pam(D: np.ndarray, k: int) -> np.ndarray:
    """PAM (build + exact swap) against a distance matrix; 0-based labels.

    ponytail: clusters must be non-empty and every point assigned to a chosen
    medoid — sane for a diverse gene table; degenerate all-identical-samples
    inputs are out of scope.
    """
    n = D.shape[0]
    # Build: greedily add the point that most reduces total distance to the
    # nearest chosen medoid.
    medoids = [int(np.argmin(D.sum(axis=1)))]
    nearest = D[:, medoids[0]].copy()
    for _ in range(k - 1):
        gain = np.maximum(nearest - D, 0).sum(axis=0)
        gain[medoids] = -1
        medoids.append(int(np.argmax(gain)))
        nearest = np.minimum(nearest, D[:, medoids[-1]])

    medoid_index = np.asarray(medoids)
    to_medoids = D[:, medoid_index].copy()
    labels = np.argmin(to_medoids, axis=1)
    while True:
        current = float(D[np.arange(n), medoid_index[labels]].sum())
        # Blank each point's own medoid so its row minimum is the best ALTERNATIVE;
        # n.inf marks points whose removal would empty their cluster (refuse swap).
        alternatives = to_medoids.copy()
        alternatives[np.arange(n), labels] = np.inf
        second_nearest = alternatives.min(axis=1)

        best_delta, best_swap = 0.0, None
        for position, medoid in enumerate(medoids):
            members = np.flatnonzero(labels == position)
            if members.size == 0 or np.isinf(second_nearest[members]).all():
                continue
            # Exact PAM pricing for moving THIS medoid to a member of its cluster.
            within = D[np.ix_(members, members)]
            deltas = (np.minimum(within, second_nearest[members][:, None]) - D[members, medoid][:, None]).sum(axis=0)
            candidate = int(np.argmin(deltas))
            if deltas[candidate] < best_delta:
                best_delta, best_swap = float(deltas[candidate]), (position, int(members[candidate]))
        if best_swap is None:
            break
        position, candidate = best_swap
        medoids[position] = candidate
        to_medoids[:, position] = D[:, candidate]
        new_labels = np.argmin(to_medoids, axis=1)
        new_cost = to_medoids[np.arange(n), new_labels].sum()
        # Commit only on real improvement: a zero-gain move would loop forever.
        if new_cost >= current - 1e-12:
            break
        labels = new_labels
    return np.argmin(D[:, np.asarray(medoids)], axis=1)


def distance_matrix(coordinates: np.ndarray) -> np.ndarray:
    """Full Euclidean distance matrix via the ||a||^2 + ||b||^2 - 2ab identity."""
    squared = (coordinates**2).sum(axis=1)
    return np.sqrt(np.maximum(squared[:, None] + squared[None, :] - 2 * coordinates @ coordinates.T, 0))


def cluster_kmedoids(coordinates: np.ndarray, k: int, random_state: int) -> np.ndarray:
    """K-medoids partition (hand-written PAM; build phase is deterministic)."""
    del random_state
    return _pam(distance_matrix(coordinates), k)


def _hierarchical_partitioner(linkage: str, metric: str) -> Callable[..., np.ndarray]:
    """Build the merge tree once, cut it at each k (the tree costs one pass)."""

    def partition(coordinates: np.ndarray, k: int, random_state: int = 0) -> np.ndarray:
        del random_state  # deterministic
        tree = scipy_linkage(coordinates, method=linkage, metric=metric)
        return fcluster(tree, k, criterion="maxclust") - 1

    return partition


#: method name -> (coordinates, k, random_state) -> 0-based labels
FIT_FUNCTIONS: dict[str, Callable[..., np.ndarray]] = {
    "kmeans": cluster_kmeans,
    "kmedoids": cluster_kmedoids,
    **{name: _hierarchical_partitioner(spec["linkage"], spec["metric"]) for name, spec in HC_METRICS.items()},
}


def fit_method(method: str, coordinates: np.ndarray, k: int, random_state: int = 42) -> np.ndarray:
    """One method's 0-based labels on the coordinate matrix; raises on unknown method."""
    if method not in FIT_FUNCTIONS:
        raise ValueError(f"Unknown clustering method: {method!r} (expected one of {METHODS})")
    return FIT_FUNCTIONS[method](coordinates, k, random_state)


# =============================================================================
# CORE LOGIC — scoring + DR-sorted labelling
# =============================================================================
def label_grid(
    coordinates: np.ndarray, methods: list[str], k_min: int, k_max: int, random_state: int
) -> pd.DataFrame:
    """All (method, k) labels in one long frame: columns method, k, gene_id, cluster (0-based)."""
    rows: list[pd.DataFrame] = []
    for method in methods:
        fit = FIT_FUNCTIONS[method]
        if method in HC_METRICS:  # one tree per method, cut at every k
            spec = HC_METRICS[method]
            tree = scipy_linkage(coordinates, method=spec["linkage"], metric=spec["metric"])
            for k in range(k_min, k_max + 1):
                partitions = fcluster(tree, k, criterion="maxclust") - 1  # 1-based -> 0-based
                rows.append(pd.DataFrame({"method": method, "k": k, "cluster": pd.Series(partitions)}))
        else:
            for k in range(k_min, k_max + 1):
                partitions = fit(coordinates, k, random_state)
                rows.append(pd.DataFrame({"method": method, "k": k, "cluster": pd.Series(partitions)}))
    return pd.concat(rows, ignore_index=True)


# =============================================================================
# CORE LOGIC — extended QC metrics (v2 k_selection lineage)
# =============================================================================
# Dunn reads only extremes (min between-gap over max diameter), so a fixed-seed
# subsample reproduces it; the full 4513² matrix would cost ~160 MB for no gain.
DUNN_MAX_POINTS = 3000
# Uniform reference draws for the gap statistic (Tibshirani's own default order).
GAP_REFERENCE_DRAWS = 50
# Bootstrap replicates for the stability ARI; k-medoids runs PAM per replicate,
# so it pays a reduced count. A flat 50 would dominate the sweep runtime.
BOOTSTRAP_REPLICATES = 50
BOOTSTRAP_REPLICATES_BY_METHOD = {"kmedoids": 20}
BOOTSTRAP_FRACTION = 0.8


def dunn_index(coordinates: np.ndarray, labels: np.ndarray, *, seed: int = 0) -> float:
    """Smallest between-cluster distance over largest within-cluster diameter."""
    if len(np.unique(labels)) < 2:
        return float("nan")
    if len(coordinates) > DUNN_MAX_POINTS:
        picked = np.random.default_rng(seed).choice(len(coordinates), size=DUNN_MAX_POINTS, replace=False)
        coordinates, labels = coordinates[picked], labels[picked]

    groups = [coordinates[labels == label] for label in np.unique(labels)]
    largest_diameter = 0.0
    for group in groups:
        if len(group) < 2:
            continue
        spread = np.sqrt(((group[:, None, :] - group[None, :, :]) ** 2).sum(axis=-1))
        largest_diameter = max(largest_diameter, float(spread.max()))
    if largest_diameter == 0:
        return float("nan")

    smallest_gap = np.inf
    for i, first in enumerate(groups):
        for second in groups[i + 1:]:
            gap = float(np.sqrt(((first[:, None, :] - second[None, :, :]) ** 2).sum(axis=-1)).min())
            smallest_gap = min(smallest_gap, gap)
    return smallest_gap / largest_diameter


def axis_shares(coordinates: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    """Fraction of the between-cluster SS carried by (DR, DL) — what the partition cuts."""
    overall = coordinates.mean(axis=0)
    between = np.zeros(coordinates.shape[1])
    for label in np.unique(labels):
        members = coordinates[labels == label]
        between += len(members) * (members.mean(axis=0) - overall) ** 2
    total = between.sum()
    if total <= 0:
        return float("nan"), float("nan")
    return float(between[0] / total), float(between[1] / total)


def score_labels(coordinates: np.ndarray, labels: np.ndarray) -> dict:
    """Full QC metric row for one partition (v2 k_selection lineage, unrounded)."""
    n_clusters = len(np.unique(labels))
    # Degenerate partitions (1 cluster, or k == n points) break every metric below;
    # NaN keeps the sweep alive and the table honest about it.
    if n_clusters < 2 or n_clusters >= len(labels):
        return {name: float("nan") for name in
                ("silhouette_score", "calinski_harabasz_score", "davies_bouldin_score", "dunn",
                 "min_cluster_silhouette", "min_cluster_share", "max_cluster_share", "n_clusters")}
    samples = silhouette_samples(coordinates, labels)
    per_cluster = np.array([samples[labels == label].mean() for label in np.unique(labels)])
    sizes = np.array([(labels == label).sum() for label in np.unique(labels)])
    return {
        "silhouette_score": float(silhouette_score(coordinates, labels)),
        "calinski_harabasz_score": float(calinski_harabasz_score(coordinates, labels)),
        "davies_bouldin_score": float(davies_bouldin_score(coordinates, labels)),
        "dunn": dunn_index(coordinates, labels),
        "min_cluster_silhouette": float(per_cluster.min()),
        "min_cluster_share": float(sizes.min() / len(labels)),
        "max_cluster_share": float(sizes.max() / len(labels)),
        "n_clusters": n_clusters,
    }


def gap_statistic(
    coordinates: np.ndarray, k_min: int, k_max: int, *, random_state: int = 42
) -> dict[int, tuple[float, float]]:
    """Tibshirani's gap at each k: observed log-inertia against a uniform null on the bounding box.

    A property of the predictors (dataset, k), not of any method — kmeans costs
    only one fit plus the null draws, so it is computed once per (dataset, k).
    """
    lows, highs = coordinates.min(axis=0), coordinates.max(axis=0)
    rng = np.random.default_rng(random_state)
    references = [rng.uniform(lows, highs, size=coordinates.shape) for _ in range(GAP_REFERENCE_DRAWS)]
    gaps: dict[int, tuple[float, float]] = {}
    for k in range(k_min, k_max + 1):
        observed = float(np.log(KMeans(n_clusters=k, n_init=10, random_state=random_state).fit(coordinates).inertia_))
        reference_logs = np.log(
            [KMeans(n_clusters=k, n_init=10, random_state=random_state).fit(ref).inertia_ for ref in references]
        )
        # Mean-of-logs is what the statistic is defined on; SE per Tibshirani 2001 eq. 4.
        gaps[k] = (float(reference_logs.mean() - observed), float(reference_logs.std(ddof=1) * np.sqrt(1 + 1 / GAP_REFERENCE_DRAWS)))
    return gaps


def partitioner(method: str, coordinates: np.ndarray, *, random_state: int = 42) -> Callable[[int], np.ndarray]:
    """A `k -> labels` closure for one method on one sample.

    Whatever the method needs before any k (a merge tree, a distance matrix) is
    paid here once — a bootstrap resample or a k sweep cuts it many times.
    """
    if method in HC_METRICS:
        spec = HC_METRICS[method]
        tree = scipy_linkage(coordinates, method=spec["linkage"], metric=spec["metric"])
        return lambda k: fcluster(tree, k, criterion="maxclust") - 1
    fit = FIT_FUNCTIONS.get(method)
    if fit is None:
        raise ValueError(f"Unknown clustering method: {method!r} (expected one of {METHODS})")
    return lambda k: fit(coordinates, k, random_state)


def bootstrap_stability(
    coordinates: np.ndarray,
    labels_by_k: dict[int, np.ndarray],
    *,
    method: str,
    random_state: int = 42,
) -> dict[int, tuple[float, float]]:
    """Mean/sd ARI of bootstrap refits vs the full-data labels, per k.

    Both loops inverted: one preparation per resample (HC tree / kmedoids
    distance matrix), then every k cut on it — hierarchical refits are the cost.
    """
    replicates = BOOTSTRAP_REPLICATES_BY_METHOD.get(method, BOOTSTRAP_REPLICATES)
    rng = np.random.default_rng(random_state)
    n = len(coordinates)
    size = int(n * BOOTSTRAP_FRACTION)
    scores: dict[int, list[float]] = {k: [] for k in labels_by_k}
    for _ in range(replicates):
        picked = rng.choice(n, size=size, replace=False)
        cut = partitioner(method, coordinates[picked], random_state=random_state)
        for k, full_labels in labels_by_k.items():
            scores[k].append(adjusted_rand_score(full_labels[picked], cut(k)))
    return {
        k: (float(np.asarray(v).mean()), float(np.asarray(v).std(ddof=1)))
        for k, v in scores.items()
    }


def sort_labels_by_dr(
    annotated: pd.DataFrame, labels: pd.Series
) -> pd.Series:
    """Renumber arbitrary labels to 1..k ascending in mean DR (depleted end = 1).

    The generic form of candidates.renumber_by_dr: the wt_cluster is the LAST id
    (highest DR, i.e. the undepleted reference sits at k), so a 9-cluster call
    with k=9 matches the final-clusters contract byte-for-byte.
    """
    k = int(labels.nunique())
    return renumber_by_dr(annotated, labels, n_clusters=k, wt_cluster=k)
