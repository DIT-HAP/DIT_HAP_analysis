#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Finalize Clusters — `direct` variant (deterministic, no human merge)
======================================================================

One of the buildable finalize variants (design doc §3.1): reads the labels produced
by cluster_one_method (already at k=final_n_clusters for a direct variant) and
deterministically renumbers clusters by DR (lowest mean DR = WT). Writes
final_clusters.tsv with the unified `cluster` column consumed by enrichment.smk +
ml.smk (design doc §4), plus a per-variant metrics.tsv (silhouette / CH / DB) for
comparing variants.

Input
-----
- annotated_data.pkl (from prepare_clustering_data)
- scaled_data.pkl (from prepare_clustering_data) — for scoring only
- labels.pkl (from cluster_one_method at k=final_n_clusters)

Output
------
- final_clusters.tsv: full annotated table + final `cluster` (1..9, WT=9);
  index = systematic ID. No raw_cluster (direct has no pre-merge labels).
- metrics.tsv: one row of silhouette / calinski_harabasz / davies_bouldin + n_clusters.

Usage
-----
    python finalize_direct_clusters.py \\
        --annotated-data results/clustering/{dataset}/_work/annotated_data.pkl \\
        --scaled-data    results/clustering/{dataset}/_work/scaled_data.pkl \\
        --labels         results/clustering/{dataset}/{variant}/_labels.pkl \\
        --output         results/clustering/final/{dataset}/{variant}/final_clusters.tsv \\
        --metrics-output results/clustering/final/{dataset}/{variant}/metrics.tsv \\
        --n-clusters 9 --wt-cluster 9

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-21
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
from workflow.src.clustering.candidates import FINAL_N_CLUSTERS, finalize_direct, score_labels


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class FinalizeDirectConfig:
    """Inputs, outputs, and params for the `direct` finalize variant."""
    annotated_data: Path
    scaled_data: Path
    labels: Path
    output: Path
    metrics_output: Path
    n_clusters: int = FINAL_N_CLUSTERS
    wt_cluster: int = 9

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.annotated_data, self.scaled_data, self.labels]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output, self.metrics_output]:
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
def run(config: FinalizeDirectConfig) -> None:
    """Renumber the direct labels by DR, write final_clusters.tsv + metrics.tsv."""
    config.validate()
    annotated = pd.read_pickle(config.annotated_data)
    scaled = pd.read_pickle(config.scaled_data)
    raw_labels = pd.read_pickle(config.labels)
    out = finalize_direct(
        annotated, raw_labels,
        n_clusters=config.n_clusters,
        wt_cluster=config.wt_cluster,
    )
    out.to_csv(config.output, sep="\t", index=True)
    logger.success(f"Wrote {len(out)} genes ({out['cluster'].notna().sum()} clustered) to {config.output}")

    # Per-variant metrics on the final clustering (label permutation-invariant scores).
    final = out.loc[scaled.index, "cluster"]
    metrics = pd.DataFrame([{"variant_type": "direct", **score_labels(scaled, final.to_numpy())}])
    metrics.to_csv(config.metrics_output, sep="\t", index=False)
    logger.success(f"Wrote metrics to {config.metrics_output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Finalize clusters — direct variant (renumber labels by DR)")
    parser.add_argument("--annotated-data", type=Path, required=True, help="Annotated data pickle (from prepare)")
    parser.add_argument("--scaled-data", type=Path, required=True, help="Scaled (DR, DL) matrix pickle (for scoring)")
    parser.add_argument("--labels", type=Path, required=True, help="Labels pickle from cluster_one_method (k=final_n_clusters)")
    parser.add_argument("--output", type=Path, required=True, help="Output final_clusters.tsv")
    parser.add_argument("--metrics-output", type=Path, required=True, help="Output per-variant metrics.tsv")
    parser.add_argument("--n-clusters", type=int, default=FINAL_N_CLUSTERS, help=f"Final cluster count (default {FINAL_N_CLUSTERS})")
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
            labels=args.labels,
            output=args.output,
            metrics_output=args.metrics_output,
            n_clusters=args.n_clusters,
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
