#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Finalize Clusters — `direct` variant (deterministic, no human merge)
======================================================================

One of the buildable finalize variants (design doc §3.1): reads the
prepare_clustering_data spine pickles, clusters the scaled (DR, DL) matrix
straight to k=9 with the chosen method, and deterministically renumbers clusters
(lowest mean DR = WT). Writes final_clusters.tsv with the unified `cluster` column
consumed by enrichment.smk + ml.smk (design doc §4).

Input
-----
- annotated_data.pkl (from prepare_clustering_data)
- scaled_data.pkl (from prepare_clustering_data)

Output
------
- final_clusters.tsv: full annotated table + final `cluster` (1..9, WT=9);
  index = systematic ID. No raw_cluster (direct has no pre-merge labels).

Usage
-----
    python finalize_direct_clusters.py \\
        --annotated-data results/clustering/candidates/{dataset}/_work/annotated_data.pkl \\
        --scaled-data    results/clustering/candidates/{dataset}/_work/scaled_data.pkl \\
        --output         results/clustering/final/{dataset}/{variant}/final_clusters.tsv \\
        --method kmeans --n-clusters 9 --random-state 42 --wt-cluster 9

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-21
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
from workflow.src.clustering.candidates import BEST_METHOD, FINAL_N_CLUSTERS, METHODS, finalize_direct


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class FinalizeDirectConfig:
    """Inputs, output, and clustering params for the `direct` finalize variant."""
    annotated_data: Path
    scaled_data: Path
    output: Path
    method: str = BEST_METHOD
    n_clusters: int = FINAL_N_CLUSTERS
    random_state: int = 42
    wt_cluster: int = 9

    def validate(self) -> None:
        """Raise ValueError on unknown method or missing input, then ensure output dir exists."""
        if self.method not in METHODS:
            raise ValueError(f"Unknown clustering method: {self.method!r} (expected one of {METHODS})")
        for path in [self.annotated_data, self.scaled_data]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


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
def run(config: FinalizeDirectConfig) -> None:
    """Cluster to k=9 with the chosen method, renumber by DR, write final_clusters.tsv."""
    config.validate()
    annotated = pd.read_pickle(config.annotated_data)
    scaled = pd.read_pickle(config.scaled_data)
    out = finalize_direct(
        annotated, scaled,
        method=config.method,
        n_clusters=config.n_clusters,
        random_state=config.random_state,
        wt_cluster=config.wt_cluster,
    )
    out.to_csv(config.output, sep="\t", index=True)
    logger.success(f"Wrote {len(out)} genes ({out['cluster'].notna().sum()} clustered) to {config.output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Finalize clusters — direct variant (k=9)")
    parser.add_argument("--annotated-data", type=Path, required=True, help="Annotated data pickle (from prepare)")
    parser.add_argument("--scaled-data", type=Path, required=True, help="Scaled (DR, DL) matrix pickle (from prepare)")
    parser.add_argument("--output", type=Path, required=True, help="Output final_clusters.tsv")
    parser.add_argument("--method", default=BEST_METHOD, choices=METHODS, help=f"Clustering method (default {BEST_METHOD})")
    parser.add_argument("--n-clusters", type=int, default=FINAL_N_CLUSTERS, help=f"Final cluster count (default {FINAL_N_CLUSTERS})")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed (default 42)")
    parser.add_argument("--wt-cluster", type=int, default=9, help="WT/background cluster id (default 9)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run finalize-direct, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = FinalizeDirectConfig(
            annotated_data=args.annotated_data,
            scaled_data=args.scaled_data,
            output=args.output,
            method=args.method,
            n_clusters=args.n_clusters,
            random_state=args.random_state,
            wt_cluster=args.wt_cluster,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
