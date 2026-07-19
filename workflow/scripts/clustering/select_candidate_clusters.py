#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Candidate Cluster Selection (aggregate)
=========================================

Final stage of the split candidate-clustering pipeline: attaches the pinned
best-method (kmeans) labels to the annotated fitting table and merges the
k-sweep + per-method metrics into one file. Reproduces the former monolithic
generate_candidate_clusters.py outputs byte-for-byte.

The best method is PINNED to kmeans (quirk #3): the source notebook's set()[0]
selection was non-deterministic but historically resolved to kmeans.

Input
-----
- annotated_data.pkl (from prepare_clustering_data)
- k_sweep_metrics.pkl (from prepare_clustering_data)
- {method}_labels.pkl for the best method (from cluster_one_method)
- {method}_metrics.pkl for all four methods (from cluster_one_method)

Output
------
- candidate_clusters.tsv: annotated table + cluster column (index = systematic ID)
- clustering_metrics.tsv: k-sweep rows + per-method comparison rows

Usage
-----
    python select_candidate_clusters.py \\
        --annotated-data .../_work/annotated_data.pkl \\
        --ksweep .../_work/k_sweep_metrics.pkl \\
        --best-labels .../_work/kmeans_labels.pkl \\
        --method-metrics .../_work/kmeans_metrics.pkl .../_work/hierarchical_agg_metrics.pkl ... \\
        --output .../candidate_clusters.tsv \\
        --metrics-output .../clustering_metrics.tsv

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
from workflow.src.clustering.candidates import BEST_METHOD, METHODS


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class SelectConfig:
    """Inputs and outputs for cluster selection + metric aggregation."""
    annotated_data: Path
    ksweep: Path
    best_labels: Path
    method_metrics: list[Path]
    output_clusters: Path
    output_metrics: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.annotated_data, self.ksweep, self.best_labels, *self.method_metrics]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_clusters, self.output_metrics]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


@logger.catch(reraise=True)
def combine_metrics(ksweep: pd.DataFrame, method_metrics: pd.DataFrame) -> pd.DataFrame:
    """Concat the k-sweep table (unrounded) and the 3dp-rounded per-method rows, byte-faithfully."""
    k_sweep_out = ksweep.assign(table="k_sweep")
    method_out = method_metrics.round(3).assign(table="method_comparison")
    return pd.concat([k_sweep_out, method_out], ignore_index=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: SelectConfig) -> None:
    """Attach best-method labels to the annotated table and write both output files."""
    config.validate()

    data_df = pd.read_pickle(config.annotated_data)
    best_labels = pd.read_pickle(config.best_labels)

    logger.info(f"Assigning cluster labels from pinned best method: {BEST_METHOD}")
    data_df["cluster"] = data_df.index.map(best_labels)
    data_df.to_csv(config.output_clusters, sep="\t")
    logger.success(f"Wrote {len(data_df)} candidate cluster labels to {config.output_clusters}")

    # Per-method metric rows, ordered canonically by METHODS regardless of CLI arg order.
    method_frames = [pd.read_pickle(p) for p in config.method_metrics]
    method_metrics = pd.concat(method_frames, ignore_index=True)
    method_metrics["method"] = pd.Categorical(method_metrics["method"], categories=METHODS, ordered=True)
    method_metrics = method_metrics.sort_values("method").reset_index(drop=True)
    method_metrics["method"] = method_metrics["method"].astype(str)

    ksweep = pd.read_pickle(config.ksweep)
    combined = combine_metrics(ksweep, method_metrics)
    combined.to_csv(config.output_metrics, sep="\t", index=False)
    logger.success(f"Wrote clustering metrics to {config.output_metrics}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Select best-method candidate clusters + aggregate metrics")
    parser.add_argument("--annotated-data", type=Path, required=True, help="Annotated data pickle (from prepare)")
    parser.add_argument("--ksweep", type=Path, required=True, help="k-sweep metrics pickle (from prepare)")
    parser.add_argument("--best-labels", type=Path, required=True, help="Best-method labels pickle (kmeans)")
    parser.add_argument("--method-metrics", type=Path, nargs="+", required=True, help="All per-method metric pickles")
    parser.add_argument("--output", type=Path, required=True, dest="output_clusters", help="Output candidate clusters tsv")
    parser.add_argument("--metrics-output", type=Path, required=True, dest="output_metrics", help="Output clustering metrics tsv")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run selection, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = SelectConfig(
            annotated_data=args.annotated_data,
            ksweep=args.ksweep,
            best_labels=args.best_labels,
            method_metrics=args.method_metrics,
            output_clusters=args.output_clusters,
            output_metrics=args.output_metrics,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
