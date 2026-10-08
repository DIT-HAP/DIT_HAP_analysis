#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Clustering grid search driver (methods x k)
============================================

Runs every clustering method over every configured k on the scaled (DR, DL)
matrix and persists three artifacts: the long label table (one row per
gene x method x k, DR-ascending cluster ids), the per-(method, k) QC metrics,
and the per-(method, k, cluster) profile (size, share, per-cluster silhouette,
means). The plots and the finalize rule read these files; this script holds
all the combinatorics.

Input
-----
- annotated_data.parquet + scaled_data.parquet (from prepare_clustering_data)

Output
------
- grid_labels.parquet: columns method, k, gene, cluster (1..k, DR-ascending)
- k_selection_metrics.tsv: one row per (method, k) with silhouette/CH/DB
- k_selection_profile.tsv: one row per (method, k, cluster)

Usage
-----
    python compute_clustering_grid.py \\
        --annotated results/clustering/{dataset}/_work/annotated_data.parquet \\
        --scaled results/clustering/{dataset}/_work/scaled_data.parquet \\
        --output-labels results/clustering/{dataset}/_work/grid_labels.parquet \\
        --output-metrics results/clustering/{dataset}/_work/k_selection_metrics.tsv \\
        --output-profile results/clustering/{dataset}/_work/k_selection_profile.tsv \\
        --methods kmeans kmedoids HC_ward HC_complete HC_average \\
        --k-min 4 --k-max 15

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-10-08
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
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from clustering.grid_search import (  # noqa: E402
    METHODS,
    axis_shares,
    bootstrap_stability,
    fit_method,
    gap_statistic,
    score_labels,
    sort_labels_by_dr,
)
from sklearn.metrics import silhouette_samples  # noqa: E402
from io_table import read_parquet, write_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class GridConfig:
    """Inputs spine parquets, output paths, and the methods x k grid."""
    annotated: Path
    scaled: Path
    output_labels: Path
    output_metrics: Path
    output_profile: Path
    methods: list[str]
    k_min: int = 4
    k_max: int = 15
    random_state: int = 42

    def validate(self) -> None:
        """Raise ValueError on missing inputs or an empty k range, then ensure output dirs exist."""
        for path in [self.annotated, self.scaled]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.k_max < self.k_min:
            raise ValueError(f"k_max ({self.k_max}) is below k_min ({self.k_min})")
        unknown = [m for m in self.methods if m not in METHODS]
        if unknown:
            raise ValueError(f"Unknown methods {unknown} (expected subset of {METHODS})")
        for out in [self.output_labels, self.output_metrics, self.output_profile]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: GridConfig) -> None:
    """Fit every (method, k), sort labels by mean DR, score, and persist the three tables."""
    config.validate()

    annotated = read_parquet(config.annotated)
    scaled = read_parquet(config.scaled)
    coordinates = scaled.to_numpy()
    logger.info(f"Scaled matrix: {scaled.shape[0]} genes x {scaled.shape[1]} features")

    # Gap statistic is a property of (dataset, k) — one sweep, broadcast to every method.
    gaps = gap_statistic(coordinates, config.k_min, config.k_max, random_state=config.random_state)

    label_rows: list[pd.DataFrame] = []
    metric_rows: list[pd.DataFrame] = []
    profile_rows: list[pd.DataFrame] = []
    for method in config.methods:
        logger.info(f"Method {method}: sweeping k={config.k_min}..{config.k_max}")
        # Raw (0-based) labels kept per k for the bootstrap refits, sorted used for output.
        raw_by_k: dict[int, np.ndarray] = {}
        sorted_by_k: dict[int, pd.Series] = {}
        for k in range(config.k_min, config.k_max + 1):
            raw_labels = pd.Series(
                fit_method(method, coordinates, k, config.random_state), index=scaled.index
            )
            # All methods get the same DR-ascending 1..k ids so cluster labels
            # are trend-comparable across methods (and k=9 matches the WT=9 contract).
            sorted_labels = sort_labels_by_dr(annotated, raw_labels)
            raw_by_k[k] = raw_labels.to_numpy()
            sorted_by_k[k] = sorted_labels

            label_rows.append(
                pd.DataFrame({"method": method, "k": k, "gene": sorted_labels.index, "cluster": sorted_labels.to_numpy()})
            )

            score = score_labels(coordinates, raw_labels.to_numpy())
            metric_rows.append(pd.DataFrame([{"method": method, "k": k, **score}]))

            # Per-cluster profile (size/share/silhouette/means) — v2-style evidence.
            sil = silhouette_samples(coordinates, raw_labels.to_numpy()) if k > 1 else np.zeros(len(raw_labels))
            prof = pd.DataFrame(
                {
                    "method": method,
                    "k": k,
                    "cluster": sorted_labels.to_numpy(),
                    "DR": annotated.loc[scaled.index, "DR"].to_numpy(),
                    "DL": annotated.loc[scaled.index, "DL"].to_numpy(),
                    "silhouette": sil,
                }
            )
            profile_rows.append(prof)

            # What the partition cuts: DR vs DL share of the between-cluster SS.
            dr_share, dl_share = axis_shares(coordinates, raw_labels.to_numpy())
            metric_rows[-1]["dr_share"] = dr_share
            metric_rows[-1]["dl_share"] = dl_share

        # Bootstrap ARI (mean/sd per k) — the stability columns of this method's rows.
        stability = bootstrap_stability(coordinates, raw_by_k, method=method, random_state=config.random_state)
        for row in metric_rows:
            if row["method"].iloc[0] == method and int(row["k"].iloc[0]) in stability:
                row["ari_mean"], row["ari_sd"] = stability[int(row["k"].iloc[0])]

    labels = pd.concat(label_rows, ignore_index=True)
    metrics = pd.concat(metric_rows, ignore_index=True)

    # Gap is the same numbers for every method at a k; broadcast into columns.
    metrics["gap"] = metrics["k"].map({k: g for k, (g, _) in gaps.items()})
    metrics["gap_sd"] = metrics["k"].map({k: s for k, (_, s) in gaps.items()})

    profile = pd.concat(profile_rows, ignore_index=True).groupby(["method", "k", "cluster"], as_index=False).agg(
        n_genes=("silhouette", "size"),
        mean_silhouette=("silhouette", "mean"),
        dr_mean=("DR", "mean"),
        dl_mean=("DL", "mean"),
    )
    profile["share"] = profile["n_genes"] / profile.groupby(["method", "k"])["n_genes"].transform("sum")

    write_parquet(labels, config.output_labels)
    metrics.to_csv(config.output_metrics, sep="\t", index=False)
    profile.to_csv(config.output_profile, sep="\t", index=False)
    logger.success(
        f"Grid done: {len(config.methods)} methods x {config.k_max - config.k_min + 1} k values -> "
        f"{len(labels)} label rows, {len(metrics)} metric rows"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Cluster the scaled DR/DL matrix for every method x k")
    parser.add_argument("--annotated", type=Path, required=True, help="annotated_data.parquet from the spine")
    parser.add_argument("--scaled", type=Path, required=True, help="scaled_data.parquet from the spine")
    parser.add_argument("--output-labels", type=Path, required=True, help="Output grid_labels.parquet")
    parser.add_argument("--output-metrics", type=Path, required=True, help="Output k_selection_metrics.tsv")
    parser.add_argument("--output-profile", type=Path, required=True, help="Output k_selection_profile.tsv")
    parser.add_argument("--methods", type=str, nargs="+", default=list(METHODS), help=f"Methods to sweep (default all of {METHODS}).")
    parser.add_argument("--k-min", type=int, default=4, help="Sweep lower bound (default 4)")
    parser.add_argument("--k-max", type=int, default=15, help="Sweep upper bound (default 15)")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed (default 42)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the grid, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = GridConfig(
            annotated=args.annotated,
            scaled=args.scaled,
            output_labels=args.output_labels,
            output_metrics=args.output_metrics,
            output_profile=args.output_profile,
            methods=args.methods,
            k_min=args.k_min,
            k_max=args.k_max,
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
