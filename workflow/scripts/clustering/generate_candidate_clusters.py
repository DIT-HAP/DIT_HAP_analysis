#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gene-Level Clustering — Candidate Labeling
============================================

Clusters genes in the 2-D depletion feature space (um = max depletion rate,
lam = lag) derived from time-resolved fitness curve fitting, producing a
candidate cluster labeling for downstream manual merge. This is the
DETERMINISTIC half of DIT_HAP_pipeline/workflow/notebooks/gene_level_clustering.ipynb
(cells 4, 6, 11, 14, 16, 18); the manual 64->9 cluster merge stays in
notebooks/clustering/finalize_gene_clusters.ipynb (design doc §5).

Input
-----
- Gene-level curve-fitting statistics (release/gene_level/fitting_results.tsv)
- Curated essentiality verification table

Output
------
- candidate_clusters.tsv: all fit columns + RevisedDeletion_essentiality + cluster
- clustering_metrics.tsv: per-method silhouette / Calinski-Harabasz / Davies-Bouldin

Usage
-----
    python generate_candidate_clusters.py \\
        --fitting-results .../release/gene_level/fitting_results.tsv \\
        --essentiality-verification-csv resources/curated/essentiality_verification.csv \\
        --output results/clustering/candidates/{dataset}/candidate_clusters.tsv \\
        --metrics-output results/clustering/candidates/{dataset}/clustering_metrics.tsv

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-16
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass, field
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
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io import read_file

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Curve-fit features considered for correlation (viz only); clustering uses um+lam.
ALL_FEATURES = ["A", "um", "lam", "t10", "t50", "t90", "t_window", "t_inflection", "y_inflection", "auc"]
SELECTED_FEATURES = ["um", "lam"]
# um above this cap is clamped; lam is divided by this divisor (byte-faithful quirk).
UM_CAP = 1.3
LAM_DIVISOR = 10
# Pinned best method — the notebook selected via set()[0] which is non-deterministic;
# it historically resolved to kmeans (hence kmeans_cluster_result.tsv). Pin it explicitly.
BEST_METHOD = "kmeans"


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class ClusteringConfig:
    """Inputs, outputs, and clustering parameters for candidate labeling."""
    fitting_results: Path
    essentiality_verification_csv: Path
    output_clusters: Path
    output_metrics: Path
    n_clusters: int = 64
    random_state: int = 42
    k_min: int = 2
    k_max: int = 20
    selected_features: list[str] = field(default_factory=lambda: list(SELECTED_FEATURES))

    def validate(self) -> None:
        """Raise ValueError if any required input is missing on disk."""
        for path in [self.fitting_results, self.essentiality_verification_csv]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.n_clusters <= 1:
            raise ValueError("n_clusters must be greater than 1.")


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )
# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_and_annotate(config: ClusteringConfig) -> pd.DataFrame:
    """Load fitting statistics and inject RevisedDeletion_essentiality at position 3."""
    data_df = read_file(config.fitting_results, index_col=[0])
    ess_df = read_file(config.essentiality_verification_csv)
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
    logger.info(f"Loaded {len(data_df)} genes from {config.fitting_results.name}")
    return data_df


@logger.catch
def scale_features(data_df: pd.DataFrame, selected_features: list[str]) -> pd.DataFrame:
    """Apply the notebook's bespoke scaling: cap um at 1.3, divide lam by 10; dropna defines the clustered set."""
    scaled_data = data_df[selected_features].dropna().copy()
    scaled_data["um"] = scaled_data["um"].apply(lambda x: x if x < UM_CAP else UM_CAP)
    scaled_data["lam"] = scaled_data["lam"].apply(lambda x: x / LAM_DIVISOR)
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


@logger.catch
def perform_clustering_analysis(data: pd.DataFrame, n_clusters: int, random_state: int) -> dict[str, np.ndarray]:
    """Apply KMeans, agglomerative-ward, divisive-complete-cityblock, and GMM; return 0-based labels per method."""
    results = {}

    logger.info("Running K-means clustering")
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10)
    results["kmeans"] = kmeans.fit_predict(data)

    logger.info("Running agglomerative hierarchical clustering")
    agg = AgglomerativeClustering(n_clusters=n_clusters, linkage="ward", metric="euclidean")
    results["hierarchical_agg"] = agg.fit_predict(data)

    logger.info("Running divisive hierarchical clustering")
    linkage_matrix = linkage(data, method="complete", metric="cityblock")
    # fcluster is 1-based; subtract 1 so all methods share 0-based labels (quirk #5).
    results["hierarchical_div"] = fcluster(linkage_matrix, n_clusters, criterion="maxclust") - 1

    logger.info("Running Gaussian Mixture Model")
    gmm = GaussianMixture(n_components=n_clusters, random_state=random_state)
    results["gmm"] = gmm.fit_predict(data)

    return results


@logger.catch
def calculate_clustering_metrics(data: pd.DataFrame, clustering_results: dict[str, np.ndarray]) -> pd.DataFrame:
    """Score each method with silhouette / Calinski-Harabasz / Davies-Bouldin, rounded to 3 dp."""
    rows = []
    for method, labels in clustering_results.items():
        rows.append(
            {
                "method": method,
                "silhouette_score": silhouette_score(data, labels),
                "calinski_harabasz_score": calinski_harabasz_score(data, labels),
                "davies_bouldin_score": davies_bouldin_score(data, labels),
                "n_clusters": len(np.unique(labels)),
            }
        )
    return pd.DataFrame(rows).round(3)


@logger.catch
def run_clustering(config: ClusteringConfig) -> pd.DataFrame:
    """Orchestrate load -> scale -> k-sweep -> multi-method clustering -> candidate labeling."""
    config.validate()
    config.output_clusters.parent.mkdir(parents=True, exist_ok=True)
    config.output_metrics.parent.mkdir(parents=True, exist_ok=True)

    data_df = load_and_annotate(config)
    scaled_data = scale_features(data_df, config.selected_features)

    k_range = range(config.k_min, config.k_max + 1)
    k_sweep = evaluate_cluster_numbers(scaled_data.values, k_range, config.random_state)

    clustering_results = perform_clustering_analysis(scaled_data, config.n_clusters, config.random_state)
    metrics_df = calculate_clustering_metrics(scaled_data, clustering_results)

    # Pin best method to kmeans (quirk #3): the notebook's set()[0] selection is non-deterministic.
    logger.info(f"Assigning cluster labels from pinned best method: {BEST_METHOD}")
    scaled_data["cluster"] = clustering_results[BEST_METHOD]
    data_df["cluster"] = data_df.index.map(scaled_data["cluster"])

    # Persist candidate labeling (index = systematic ID) and both metric tables.
    data_df.to_csv(config.output_clusters, sep="\t")
    logger.success(f"Wrote {len(data_df)} candidate cluster labels to {config.output_clusters}")

    # Merge the k-sweep and per-method metrics into one metrics file (k-sweep rows + method rows).
    k_sweep_out = k_sweep.assign(table="k_sweep")
    method_out = metrics_df.assign(table="method_comparison")
    combined = pd.concat([k_sweep_out, method_out], ignore_index=True)
    combined.to_csv(config.output_metrics, sep="\t", index=False)
    logger.success(f"Wrote clustering metrics to {config.output_metrics}")

    return data_df
# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Generate candidate gene-level clusters (deterministic)")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Gene-level curve-fitting statistics tsv")
    parser.add_argument("--essentiality-verification-csv", type=Path, required=True, help="Curated essentiality verification csv")
    parser.add_argument("--output", type=Path, required=True, dest="output_clusters", help="Output candidate clusters tsv")
    parser.add_argument("--metrics-output", type=Path, required=True, dest="output_metrics", help="Output clustering metrics tsv")
    parser.add_argument("--n-clusters", type=int, default=64, help="Number of candidate clusters (default 64)")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed (default 42)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run clustering, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = ClusteringConfig(
            fitting_results=args.fitting_results,
            essentiality_verification_csv=args.essentiality_verification_csv,
            output_clusters=args.output_clusters,
            output_metrics=args.output_metrics,
            n_clusters=args.n_clusters,
            random_state=args.random_state,
        )
        run_clustering(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
