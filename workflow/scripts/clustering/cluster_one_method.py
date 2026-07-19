#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Single-Method Clustering
==========================

Runs ONE clustering method (kmeans / hierarchical_agg / hierarchical_div / gmm)
over the scaled (DR, DL) matrix from prepare_clustering_data, and emits that
method's 0-based labels + metric row. Fanned out per method by clustering.smk's
`method` wildcard (design doc §5) — the analogue of ml.smk's target x mode split.

All four methods are deterministic given the scaled matrix + seed, so running
each in its own process reproduces the monolith's results exactly.

Input
-----
- scaled_data.pkl: the scaled (DR, DL) matrix (from prepare_clustering_data)

Output
------
- {method}_labels.pkl: pd.Series of 0-based cluster labels, indexed by systematic ID
- {method}_metrics.pkl: one-row DataFrame (method + silhouette/CH/DB + n_clusters, unrounded)

Usage
-----
    python cluster_one_method.py \\
        --method kmeans \\
        --scaled-data results/clustering/candidates/{dataset}/_work/scaled_data.pkl \\
        --output-labels results/clustering/candidates/{dataset}/_work/kmeans_labels.pkl \\
        --output-metrics results/clustering/candidates/{dataset}/_work/kmeans_metrics.pkl

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-17
Version:  1.0.0
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
from workflow.src.clustering.candidates import METHODS, cluster_one_method, score_labels


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class MethodConfig:
    """Inputs, outputs, and parameters for a single clustering method."""
    method: str
    scaled_data: Path
    output_labels: Path
    output_metrics: Path
    n_clusters: int = 64
    random_state: int = 42

    def validate(self) -> None:
        """Raise ValueError on an unknown method, missing input, or bad n_clusters; ensure output dirs exist."""
        if self.method not in METHODS:
            raise ValueError(f"Unknown clustering method: {self.method!r} (expected one of {METHODS})")
        if not self.scaled_data.exists():
            raise ValueError(f"Required input not found: {self.scaled_data}")
        if self.n_clusters <= 1:
            raise ValueError("n_clusters must be greater than 1.")
        for out in [self.output_labels, self.output_metrics]:
            out.parent.mkdir(parents=True, exist_ok=True)


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
    """Run one clustering method on the scaled matrix and persist its labels + metric row."""
    config.validate()
    scaled_data = pd.read_pickle(config.scaled_data)

    labels = cluster_one_method(config.method, scaled_data, config.n_clusters, config.random_state)
    # Index labels by systematic ID so select_candidate_clusters can map them back safely.
    label_series = pd.Series(labels, index=scaled_data.index, name="cluster")
    label_series.to_pickle(config.output_labels)

    metrics_row = pd.DataFrame([{"method": config.method, **score_labels(scaled_data, labels)}])
    metrics_row.to_pickle(config.output_metrics)
    logger.success(f"[{config.method}] {len(label_series)} genes labeled -> {config.output_labels}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Run one gene-level clustering method (deterministic)")
    parser.add_argument("--method", required=True, choices=METHODS, help="Clustering method to run")
    parser.add_argument("--scaled-data", type=Path, required=True, help="Scaled (DR, DL) matrix pickle")
    parser.add_argument("--output-labels", type=Path, required=True, help="Output labels pickle (pd.Series)")
    parser.add_argument("--output-metrics", type=Path, required=True, help="Output one-row metrics pickle")
    parser.add_argument("--n-clusters", type=int, default=64, help="Number of candidate clusters (default 64)")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed (default 42)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run one method, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = MethodConfig(
            method=args.method,
            scaled_data=args.scaled_data,
            output_labels=args.output_labels,
            output_metrics=args.output_metrics,
            n_clusters=args.n_clusters,
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
