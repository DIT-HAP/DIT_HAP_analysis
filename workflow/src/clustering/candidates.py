"""
Gene-Level Clustering — Candidate Helpers
===========================================

Shared logic for the split candidate-clustering pipeline: constants, the
load/annotate/scale/k-sweep preprocessing, per-method clustering dispatch, and
metric scoring. Used by the three driver scripts
(prepare_clustering_data / cluster_one_method / select_candidate_clusters) so
each clustering method can be its own Snakemake job while the maths lives in one
importable place.

The four methods all consume the SAME scaled (DR, DL) matrix and params; they
differ only in the estimator — hence one parameterized `cluster_one_method`
dispatched by name rather than four near-identical functions.

Input
-----
- Gene-level curve-fitting statistics + curated essentiality verification (prepare)
- The scaled (DR, DL) matrix pickle (per-method + select)

Output
------
- Helpers returning DataFrames / label arrays; drivers persist them.

Usage
-----
    from workflow.src.clustering.candidates import scale_features, cluster_one_method
    labels = cluster_one_method("kmeans", scaled_data, n_clusters=64, random_state=42)

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-17
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.metrics import (
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_score,
)
from sklearn.mixture import GaussianMixture

# 4. Local Imports
from workflow.src.io import read_file

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Curve-fit features considered for correlation (viz only); clustering uses DR+DL.
ALL_FEATURES = ["A", "DR", "DL", "t10", "t50", "t90", "t_window", "t_inflection", "y_inflection", "auc"]
SELECTED_FEATURES = ["DR", "DL"]
# DR above this cap is clamped; DL is divided by this divisor (byte-faithful quirk).
DR_CAP = 1.3
DL_DIVISOR = 10
# The four clustering methods, in the fixed order the metrics table reports them.
METHODS = ["kmeans", "hierarchical_agg", "hierarchical_div", "gmm"]
# Pinned best method — the notebook selected via set()[0] which is non-deterministic;
# it historically resolved to kmeans (hence kmeans_cluster_result.tsv). Pin it explicitly.
BEST_METHOD = "kmeans"


# =============================================================================
# PREPROCESSING (the "spine")
# =============================================================================
# Legacy -> current metric column names. Commit 573aafd renamed the clustering
# feature columns um->DR (max depletion rate) and lam->DL (lag) in code, but some
# upstream fitting_results.tsv exports still carry the old um/lam headers. Rename
# on load so the pipeline speaks DR/DL consistently; a no-op once upstream ships
# DR/DL directly (the two are the same metrics, only renamed).
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}


@logger.catch
def load_and_annotate(fitting_results: Path, essentiality_verification_csv: Path) -> pd.DataFrame:
    """Load fitting statistics, normalize legacy um/lam -> DR/DL, and inject RevisedDeletion_essentiality at position 3."""
    data_df = read_file(fitting_results, index_col=[0])
    rename = {old: new for old, new in _LEGACY_METRIC_RENAME.items() if old in data_df.columns and new not in data_df.columns}
    if rename:
        logger.info(f"Normalizing legacy metric columns: {rename}")
        data_df = data_df.rename(columns=rename)
    ess_df = read_file(essentiality_verification_csv)
    verification_essentiality = ess_df.set_index("systematic_id")["verification_essentiality"].to_dict()

    # Byte-faithful: verified value if present, else fall back to DeletionLibrary_essentiality.
    data_df.insert(
        3,
        "RevisedDeletion_essentiality",
        data_df.apply(
            lambda row: verification_essentiality[row.name]
            if row.name in verification_essentiality
            else row["DeletionLibrary_essentiality"],
            axis=1,
        ),
    )
    logger.info(f"Loaded {len(data_df)} genes from {fitting_results.name}")
    return data_df


@logger.catch
def scale_features(
    data_df: pd.DataFrame, selected_features: list[str], dr_cap: float = DR_CAP, dl_divisor: float = DL_DIVISOR
) -> pd.DataFrame:
    """Apply the notebook's bespoke scaling: cap DR at dr_cap, divide DL by dl_divisor; dropna defines the clustered set."""
    scaled_data = data_df[selected_features].dropna().copy()
    scaled_data["DR"] = scaled_data["DR"].apply(lambda x: x if x < dr_cap else dr_cap)
    scaled_data["DL"] = scaled_data["DL"].apply(lambda x: x / dl_divisor)
    logger.info(f"Scaled feature matrix: {scaled_data.shape[0]} genes x {scaled_data.shape[1]} features")
    return scaled_data


@logger.catch
def evaluate_cluster_numbers(data: np.ndarray, k_range: range, random_state: int) -> pd.DataFrame:
    """Sweep k over k_range with KMeans, recording inertia + silhouette/CH/DB metrics."""
    rows = []
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=random_state, n_init=10)
        labels = kmeans.fit_predict(data)
        rows.append(
            {
                "k": k,
                "inertia": kmeans.inertia_,
                "silhouette": silhouette_score(data, labels),
                "calinski_harabasz": calinski_harabasz_score(data, labels),
                "davies_bouldin": davies_bouldin_score(data, labels),
            }
        )
    return pd.DataFrame(rows)


# =============================================================================
# PER-METHOD CLUSTERING
# =============================================================================
@logger.catch(reraise=True)
def cluster_one_method(method: str, data: pd.DataFrame, n_clusters: int, random_state: int) -> np.ndarray:
    """Run one clustering method by name, returning 0-based labels aligned to `data`'s rows."""
    if method == "kmeans":
        logger.info("Running K-means clustering")
        return KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10).fit_predict(data)
    if method == "hierarchical_agg":
        logger.info("Running agglomerative hierarchical clustering")
        return AgglomerativeClustering(n_clusters=n_clusters, linkage="ward", metric="euclidean").fit_predict(data)
    if method == "hierarchical_div":
        logger.info("Running divisive hierarchical clustering")
        linkage_matrix = linkage(data, method="complete", metric="cityblock")
        # fcluster is 1-based; subtract 1 so all methods share 0-based labels (quirk #5).
        return fcluster(linkage_matrix, n_clusters, criterion="maxclust") - 1
    if method == "gmm":
        logger.info("Running Gaussian Mixture Model")
        return GaussianMixture(n_components=n_clusters, random_state=random_state).fit_predict(data)
    raise ValueError(f"Unknown clustering method: {method!r} (expected one of {METHODS})")


@logger.catch
def score_labels(data: pd.DataFrame, labels: np.ndarray) -> dict:
    """Score one labeling with silhouette / Calinski-Harabasz / Davies-Bouldin + cluster count (unrounded)."""
    return {
        "silhouette_score": silhouette_score(data, labels),
        "calinski_harabasz_score": calinski_harabasz_score(data, labels),
        "davies_bouldin_score": davies_bouldin_score(data, labels),
        "n_clusters": len(np.unique(labels)),
    }


# =============================================================================
# FINALIZE VARIANTS (deterministic 64/scaled -> 9, no human merge)
# =============================================================================
# Number of final clusters for the automatic finalize paths (design doc §2).
FINAL_N_CLUSTERS = 9


@logger.catch(reraise=True)
def renumber_by_dr(
    annotated: pd.DataFrame,
    raw_labels: pd.Series,
    n_clusters: int = FINAL_N_CLUSTERS,
    wt_cluster: int = 9,
) -> pd.Series:
    """Deterministically renumber arbitrary group labels to 1..n_clusters by mean DR.

    Every finalize variant (direct / auto_merge / grid / manual) ends here so the
    numbering is identical and cross-variant comparable (design doc §2): the group
    with the lowest mean DR becomes wt_cluster (WT), the rest ascend in mean-DR
    order (ties broken by mean DL then the group's own label). `raw_labels` is a
    Series of group labels indexed by gene id; returns the final 1..n_clusters
    labels on that same index. Raises if the group count != n_clusters.
    """
    stats = (
        annotated.loc[raw_labels.index, ["DR", "DL"]]
        .assign(_raw=raw_labels.to_numpy())
        .groupby("_raw")
        .agg(mean_dr=("DR", "mean"), mean_dl=("DL", "mean"))
        .reset_index()
        .sort_values(["mean_dr", "mean_dl", "_raw"], kind="stable")
        .reset_index(drop=True)
    )
    # Ascending DR -> the lowest-DR group becomes wt_cluster; the remaining ids
    # 1..n_clusters (excluding wt) fill the other ranks in ascending-DR order.
    remaining = [i for i in range(1, n_clusters + 1) if i != wt_cluster]
    ordered_ids = [wt_cluster] + remaining
    # Guard against silent mis-numbering: dict(zip(...)) would quietly truncate if
    # counts diverge (too few groups, or wt_cluster out of the 1..n_clusters range).
    if len(stats) != len(ordered_ids):
        raise ValueError(
            f"Expected {len(ordered_ids)} groups but got {len(stats)}; "
            f"cannot renumber deterministically (check n_clusters / merge / grid cuts)."
        )
    # stats is already sorted ascending-DR and re-indexed, so row order == rank order:
    # rank 0 (lowest DR) -> ordered_ids[0] (wt_cluster), rank 1 -> id 1, ...
    raw_to_final = dict(zip(stats["_raw"], ordered_ids))
    return raw_labels.map(raw_to_final)


@logger.catch(reraise=True)
def finalize_direct(
    annotated: pd.DataFrame,
    raw_labels: pd.Series,
    n_clusters: int = FINAL_N_CLUSTERS,
    wt_cluster: int = 9,
) -> pd.DataFrame:
    """`direct` variant: `raw_labels` are already at n_clusters (from cluster_one_method
    at k=final_n_clusters); just renumber by DR. Returns the annotated table with a
    final `cluster` column (NaN for genes not in the clustered set). See design doc §3.1.
    """
    out = annotated.copy()
    out["cluster"] = renumber_by_dr(annotated, raw_labels, n_clusters, wt_cluster)
    logger.info(f"finalize_direct: {len(raw_labels)} genes -> {n_clusters} clusters (WT={wt_cluster})")
    return out


@logger.catch(reraise=True)
def finalize_auto_merge(
    annotated: pd.DataFrame,
    scaled: pd.DataFrame,
    raw_labels: pd.Series,
    n_clusters: int = FINAL_N_CLUSTERS,
    wt_cluster: int = 9,
) -> pd.DataFrame:
    """`auto_merge` variant: ward-merge the transitional-cluster CENTROIDS down to
    n_clusters groups, then renumber by DR. `raw_labels` are the per-gene transitional
    labels (from cluster_one_method at k=n_intermediate). Centroids are the unweighted
    per-cluster means in scaled (DR, DL) space. Returns the annotated table with final
    `cluster` + a preserved `raw_cluster` (the transitional labels). See design doc §3.2.
    """
    # One centroid per transitional cluster, ordered by label for determinism.
    centroids = scaled.loc[raw_labels.index].groupby(raw_labels.to_numpy()).mean().sort_index()
    if len(centroids) < n_clusters:
        raise ValueError(
            f"auto_merge needs >= {n_clusters} transitional clusters to merge, got {len(centroids)}."
        )
    merge_labels = AgglomerativeClustering(n_clusters=n_clusters, linkage="ward", metric="euclidean").fit_predict(
        centroids.to_numpy()
    )
    candidate_to_group = dict(zip(centroids.index, merge_labels))
    merged = raw_labels.map(candidate_to_group)

    out = annotated.copy()
    out["cluster"] = renumber_by_dr(annotated, merged, n_clusters, wt_cluster)
    out["raw_cluster"] = raw_labels
    logger.info(
        f"finalize_auto_merge: {len(centroids)} transitional clusters -> {n_clusters} (WT={wt_cluster})"
    )
    return out


@logger.catch(reraise=True)
def finalize_grid(
    annotated: pd.DataFrame,
    scaled: pd.DataFrame,
    dr_cuts: list[float],
    dl_cuts: list[float],
    n_clusters: int = FINAL_N_CLUSTERS,
    wt_cluster: int = 9,
) -> pd.DataFrame:
    """`grid` variant: split the scaled DR/DL axes at dr_cuts / dl_cuts into a
    (len(dr_cuts)+1) x (len(dl_cuts)+1) rectangular grid, assign each gene to its
    cell, then renumber by DR. Cell count must equal n_clusters. Cuts are thresholds
    in SCALED space (DR already capped, DL already divided). See design doc §3.3.
    """
    n_cells = (len(dr_cuts) + 1) * (len(dl_cuts) + 1)
    if n_cells != n_clusters:
        raise ValueError(
            f"grid produces {len(dr_cuts) + 1}x{len(dl_cuts) + 1} = {n_cells} cells "
            f"but final_n_clusters is {n_clusters}; adjust dr_cuts / dl_cuts."
        )
    dr_bin = np.digitize(scaled["DR"].to_numpy(), sorted(dr_cuts))
    dl_bin = np.digitize(scaled["DL"].to_numpy(), sorted(dl_cuts))
    cell = pd.Series(dr_bin * (len(dl_cuts) + 1) + dl_bin, index=scaled.index, name="_cell")

    out = annotated.copy()
    out["cluster"] = renumber_by_dr(annotated, cell, n_clusters, wt_cluster)
    logger.info(f"finalize_grid: {n_cells}-cell grid -> {n_clusters} clusters (WT={wt_cluster})")
    return out
