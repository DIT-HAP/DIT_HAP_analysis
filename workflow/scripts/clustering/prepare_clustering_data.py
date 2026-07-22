#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Clustering Preprocessing (the spine)
======================================

First stage of the split candidate-clustering pipeline: loads the gene-level
fitting statistics, injects RevisedDeletion_essentiality, applies the bespoke
(DR, DL) scaling, and runs the KMeans k-sweep. Emits the annotated table + the
scaled matrix as pickles for the per-method clustering jobs and the final
select step (design doc §5).

Input
-----
- Gene-level curve-fitting statistics (release/gene_level/fitting_results.tsv)
- Curated essentiality verification table

Output
------
- annotated_data.parquet: full fitting table + RevisedDeletion_essentiality (index = systematic ID)
- scaled_data.parquet: the scaled (DR, DL) matrix that defines the clustered gene set
- k_sweep_metrics.parquet: KMeans k-sweep (inertia + silhouette/CH/DB per k)

Usage
-----
    python prepare_clustering_data.py \\
        --fitting-results .../release/gene_level/fitting_results.tsv \\
        --essentiality-verification-csv resources/curated/essentiality_verification.csv \\
        --output-annotated results/clustering/{dataset}/_work/annotated_data.parquet \\
        --output-scaled results/clustering/{dataset}/_work/scaled_data.parquet \\
        --output-ksweep results/clustering/{dataset}/_work/k_sweep_metrics.parquet

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
from dataclasses import dataclass, field
from pathlib import Path

# 2. Third-party Imports
from loguru import logger

# 3. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.clustering.candidates import (
    DL_DIVISOR,
    DR_CAP,
    SELECTED_FEATURES,
    evaluate_cluster_numbers,
    load_and_annotate,
    scale_features,
)


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PrepareConfig:
    """Inputs, outputs, and k-sweep parameters for clustering preprocessing."""
    fitting_results: Path
    essentiality_verification_csv: Path
    output_annotated: Path
    output_scaled: Path
    output_ksweep: Path
    random_state: int = 42
    k_min: int = 2
    k_max: int = 20
    dr_cap: float = DR_CAP
    dl_divisor: float = DL_DIVISOR
    selected_features: list[str] = field(default_factory=lambda: list(SELECTED_FEATURES))

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.fitting_results, self.essentiality_verification_csv]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_annotated, self.output_scaled, self.output_ksweep]:
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
def run(config: PrepareConfig) -> None:
    """Load -> annotate -> scale -> k-sweep, then persist the three pickles."""
    config.validate()

    data_df = load_and_annotate(config.fitting_results, config.essentiality_verification_csv)
    scaled_data = scale_features(data_df, config.selected_features, dr_cap=config.dr_cap, dl_divisor=config.dl_divisor)

    k_range = range(config.k_min, config.k_max + 1)
    k_sweep = evaluate_cluster_numbers(scaled_data.values, k_range, config.random_state)

    write_parquet(data_df, config.output_annotated)
    write_parquet(scaled_data, config.output_scaled)
    write_parquet(k_sweep, config.output_ksweep)
    logger.success(
        f"Prepared clustering data: {len(data_df)} genes, {len(scaled_data)} clustered, "
        f"k-sweep {config.k_min}..{config.k_max}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare clustering data (spine: load/annotate/scale/k-sweep)")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Gene-level curve-fitting statistics tsv")
    parser.add_argument("--essentiality-verification-csv", type=Path, required=True, help="Curated essentiality verification csv")
    parser.add_argument("--output-annotated", type=Path, required=True, help="Output annotated data pickle")
    parser.add_argument("--output-scaled", type=Path, required=True, help="Output scaled (DR, DL) matrix pickle")
    parser.add_argument("--output-ksweep", type=Path, required=True, help="Output k-sweep metrics pickle")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed (default 42)")
    parser.add_argument("--k-min", type=int, default=2, help="k-sweep lower bound (default 2)")
    parser.add_argument("--k-max", type=int, default=20, help="k-sweep upper bound (default 20)")
    parser.add_argument("--dr-cap", type=float, default=DR_CAP, help=f"DR clamp ceiling (default {DR_CAP})")
    parser.add_argument("--dl-divisor", type=float, default=DL_DIVISOR, help=f"DL scaling divisor (default {DL_DIVISOR})")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run preprocessing, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PrepareConfig(
            fitting_results=args.fitting_results,
            essentiality_verification_csv=args.essentiality_verification_csv,
            output_annotated=args.output_annotated,
            output_scaled=args.output_scaled,
            output_ksweep=args.output_ksweep,
            random_state=args.random_state,
            k_min=args.k_min,
            k_max=args.k_max,
            dr_cap=args.dr_cap,
            dl_divisor=args.dl_divisor,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
