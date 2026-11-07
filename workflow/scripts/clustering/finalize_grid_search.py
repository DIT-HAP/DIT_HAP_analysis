#!/usr/bin/env python3

"""
Finalize Grid Search Clusters — the selected (method, k) partition
==================================================================

Reads the grid-search label table and emits ONE final clusters table: the
configured method at the configured k, renumbered 1..k ascending in mean DR
(cluster 1 = most depleted, k = the WT/undepleted reference), with genes the
space does not cover left as NaN — the same `final_clusters.tsv` contract the
variants build, so enrichment / ml downstream read it unchanged. Also writes
the selected partition's QC metrics beside it.

Input
-----
- grid_labels.parquet (from compute_clustering_grid)
- annotated_data.parquet (keep the full annotated table + Raw DR/DL)

Output
------
- final_clusters.tsv: annotated table + `cluster` column (1..k, WT = k)
- metrics.tsv: the selected (method, k) silhouette/CH/DB row + selection

Usage
-----
    python finalize_grid_search.py \\
        --labels results/3b_clustering/{dataset}/_work/grid_labels.parquet \\
        --metrics results/3b_clustering/{dataset}/_work/k_selection_metrics.tsv \\
        --annotated results/3b_clustering/{dataset}/_work/annotated_data.parquet \\
        --output results/3b_clustering/{dataset}/{variant}/final_clusters.tsv \\
        --metrics-output results/3b_clustering/{dataset}/{variant}/metrics.tsv \\
        --method kmeans --k 9

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-10-08
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

# 3. Third-party Imports
from loguru import logger

# 2. Data Processing Imports
import pandas as pd

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from io_table import read_parquet
from logging_setup import setup_logger


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class FinalizeConfig:
    """Grid label table + annotated table in; the selected (method, k) out."""

    labels: Path
    metrics: Path
    annotated: Path
    output: Path
    metrics_output: Path
    method: str
    k: int

    def validate(self) -> None:
        """Raise ValueError on missing inputs, then ensure output dirs exist."""
        for path in [self.labels, self.metrics, self.annotated]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)
        self.metrics_output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: FinalizeConfig) -> None:
    """Pick (method, k) labels from the grid, rename to `cluster`, write final_clusters.tsv +
    metrics.tsv."""
    config.validate()

    labels = pd.read_parquet(config.labels)
    selected = labels[(labels["method"] == config.method) & (labels["k"] == config.k)]
    if selected.empty:
        available = sorted(map(tuple, labels[["method", "k"]].drop_duplicates().to_numpy()))
        raise ValueError(
            f"Grid has no labels for method={config.method!r} k={config.k}; available: {available}"
        )

    annotated = read_parquet(config.annotated)
    annotated["cluster"] = annotated.index.map(selected.set_index("gene")["cluster"])

    config.output.parent.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(config.output, sep="\t", index=True)

    metrics = pd.read_csv(config.metrics, sep="\t")
    row = metrics[(metrics["method"] == config.method) & (metrics["k"] == config.k)].copy()
    row.insert(0, "selection", f"{config.method}_k{config.k}")
    row.to_csv(config.metrics_output, sep="\t", index=False)

    n_clustered = annotated["cluster"].notna().sum()
    logger.success(
        f"finalize_grid_search: {config.method} k={config.k} -> {n_clustered} genes, "
        f"{annotated['cluster'].dropna().nunique()} clusters ({config.output})"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Emit the selected (method, k) final_clusters.tsv")
    parser.add_argument(
        "--labels",
        type=Path,
        required=True,
        help="grid_labels.parquet from compute_clustering_grid",
    )
    parser.add_argument(
        "--metrics",
        type=Path,
        required=True,
        help="k_selection_metrics.tsv from compute_clustering_grid",
    )
    parser.add_argument(
        "--annotated", type=Path, required=True, help="annotated_data.parquet from the spine"
    )
    parser.add_argument("--output", type=Path, required=True, help="Output final_clusters.tsv")
    parser.add_argument("--metrics-output", type=Path, required=True, help="Output metrics.tsv")
    parser.add_argument(
        "--method", type=str, required=True, help="Selected method name (e.g. kmeans, HC_ward)"
    )
    parser.add_argument("--k", type=int, required=True, help="Selected cluster count")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging"
    )
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, write outputs, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = FinalizeConfig(
            labels=args.labels,
            metrics=args.metrics,
            annotated=args.annotated,
            output=args.output,
            metrics_output=args.metrics_output,
            method=args.method,
            k=args.k,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
