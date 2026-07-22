#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Per-Variant Clustering (label generation)
===========================================

Runs ONE clustering method (kmeans / hierarchical_agg / hierarchical_div / gmm)
over the scaled (DR, DL) matrix from prepare_clustering_data and emits that
variant's 0-based labels. Clusters to `n_intermediate` when the variant declares
one (merge variants) else straight to `final_n_clusters` (direct variants). Fanned
out per variant by clustering.smk (design doc 2026-07-21-clustering-finalize-variants).

For merge variants the labels are the transitional grouping a later merge rule
reduces to final_n_clusters; for direct variants they are already at final_n_clusters
and only need DR-renumbering downstream.

Input
-----
- scaled_data.parquet: the scaled (DR, DL) matrix (from prepare_clustering_data)

Output
------
- labels.parquet: pd.Series of 0-based cluster labels, indexed by systematic ID

Usage
-----
    python cluster_one_method.py \\
        --method kmeans \\
        --scaled-data results/clustering/{dataset}/_work/scaled_data.parquet \\
        --output-labels results/clustering/{dataset}/{variant}/_labels.parquet \\
        --final-n-clusters 9 --n-intermediate 64 --random-state 42

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-17
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io import read_parquet, write_parquet
from workflow.src.clustering.candidates import FINAL_N_CLUSTERS, METHODS, cluster_one_method


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class MethodConfig:
    """Inputs, output, and parameters for one variant's clustering.

    Clusters to `n_intermediate` when set (merge variants), else `final_n_clusters`
    (direct variants).
    """
    method: str
    scaled_data: Path
    output_labels: Path
    final_n_clusters: int = FINAL_N_CLUSTERS
    n_intermediate: int | None = None
    random_state: int = 42

    @property
    def n_clusters(self) -> int:
        """Effective k: the transitional count if given, else the final count."""
        return self.n_intermediate if self.n_intermediate is not None else self.final_n_clusters

    def validate(self) -> None:
        """Raise ValueError on an unknown method, missing input, or bad k; ensure output dir exists."""
        if self.method not in METHODS:
            raise ValueError(f"Unknown clustering method: {self.method!r} (expected one of {METHODS})")
        if not self.scaled_data.exists():
            raise ValueError(f"Required input not found: {self.scaled_data}")
        if self.n_clusters <= 1:
            raise ValueError("n_clusters must be greater than 1.")
        if self.n_intermediate is not None and self.n_intermediate < self.final_n_clusters:
            raise ValueError(
                f"n_intermediate ({self.n_intermediate}) must be >= final_n_clusters ({self.final_n_clusters})"
            )
        self.output_labels.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: MethodConfig) -> None:
    """Cluster the scaled matrix to the effective k and persist the labels."""
    config.validate()
    scaled_data = read_parquet(config.scaled_data)

    labels = cluster_one_method(config.method, scaled_data, config.n_clusters, config.random_state)
    # Index labels by systematic ID so the finalize/merge rules can map them back safely.
    label_series = pd.Series(labels, index=scaled_data.index, name="cluster")
    write_parquet(label_series, config.output_labels)
    logger.success(f"[{config.method}] {len(label_series)} genes labeled at k={config.n_clusters} -> {config.output_labels}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Cluster one variant's method (labels only)")
    parser.add_argument("--method", required=True, choices=METHODS, help="Clustering method to run")
    parser.add_argument("--scaled-data", type=Path, required=True, help="Scaled (DR, DL) matrix pickle")
    parser.add_argument("--output-labels", type=Path, required=True, help="Output labels parquet (pd.Series)")
    parser.add_argument("--final-n-clusters", type=int, default=FINAL_N_CLUSTERS, help=f"Final cluster count (default {FINAL_N_CLUSTERS})")
    parser.add_argument("--n-intermediate", type=int, default=None, help="Transitional cluster count for merge variants (default: none -> cluster to final)")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed (default 42)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run one variant's clustering, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = MethodConfig(
            method=args.method,
            scaled_data=args.scaled_data,
            output_labels=args.output_labels,
            final_n_clusters=args.final_n_clusters,
            n_intermediate=args.n_intermediate,
            random_state=args.random_state,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
