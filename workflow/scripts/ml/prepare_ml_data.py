#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ML Modeling-Data Preparation (the spine)
==========================================

Merges the per-gene feature matrix with curve-fit targets + cluster labels and
applies the DR > threshold filter ONCE, emitting a modeling_data parquet shared
by the four target x mode AutoML jobs (which previously each re-merged the same
data). Deterministic and target/mode-independent.

Input
-----
- Per-gene feature matrix (results/features/{version}/pombe_coding_gene_protein_features.tsv)
- Curated final_clusters.tsv (Systematic ID, A, DR, DL, cluster)

Output
------
- modeling_data.parquet: merged, DR-filtered modeling table

Usage
-----
    python prepare_ml_data.py \\
        --feature-matrix results/features/2025-10-01/pombe_coding_gene_protein_features.tsv \\
        --final-clusters resources/curated/final_clusters.tsv \\
        --output results/ml/models/{dataset}/{version}/_work/modeling_data.parquet

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

# 2. Third-party Imports
from loguru import logger

# 3. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io import read_parquet, write_parquet
from workflow.src.ml.data import DR_FILTER, load_modeling_data


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PrepareConfig:
    """Inputs, output, and the DR filter for modeling-data preparation."""
    feature_matrix: Path
    final_clusters: Path
    output: Path
    dr_filter: float = DR_FILTER

    def validate(self) -> None:
        """Raise ValueError if a required input is missing, then ensure the output dir exists."""
        for path in [self.feature_matrix, self.final_clusters]:
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
def run(config: PrepareConfig) -> None:
    """Merge + filter the modeling data once and parquet it for the AutoML jobs."""
    config.validate()
    data = load_modeling_data(config.feature_matrix, config.final_clusters, config.dr_filter)
    write_parquet(data, config.output)
    logger.success(f"Wrote modeling data {data.shape} -> {config.output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare shared ML modeling data (merge + DR filter)")
    parser.add_argument("--feature-matrix", type=Path, required=True, help="Per-gene feature matrix tsv")
    parser.add_argument("--final-clusters", type=Path, required=True, help="Curated final_clusters.tsv")
    parser.add_argument("--output", type=Path, required=True, help="Output modeling_data pickle")
    parser.add_argument("--dr-filter", type=float, default=DR_FILTER, help="DR > filter threshold (default 0.3)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, prepare modeling data, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PrepareConfig(
            feature_matrix=args.feature_matrix,
            final_clusters=args.final_clusters,
            output=args.output,
            dr_filter=args.dr_filter,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
