#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare Fitness Table
======================

Stage 1 of the comparison split: load the curated DIT-HAP gene clusters, the
PomBase-derived protein-features table, and the gRNA (HD data) fitted
parameters, then merge them into the fitness_table.parquet intermediate
consumed by the stats / figures rules.

Input
-----
- final_clusters.tsv (Systematic ID + the DIT-HAP fitness metric; DR is the
  current name, legacy releases ship it as `um`, normalized on load) from the
  clustering finalize-variant stage, sourced via final_clusters_path(dataset,
  selected_variant). Only Systematic ID + DR are read here.
- gRNA HD-data fitted parameters TSV (Systematic ID + `um` gRNA fitness metric).
- pombe_coding_gene_protein_features.tsv (gene_systematic_id + the other
  large-scale study columns: Barseq_from_dulab/koch, integration density, ipkm,
  uipkm, colony_size_Malecki2016, Max Growth Rate, Colony Formation).

Output
------
- fitness_table.parquet: protein-features spine left-joined with the DIT-HAP
  (um_DIT_HAP) and gRNA (um_gRNA) fitness metrics, integration-density columns
  clipped at --clip-upper.

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-22
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
from workflow.src.io import write_parquet  # noqa: E402
from workflow.src.comparison.core import (  # noqa: E402
    CLIP_UPPER,
    build_fitness_table,
    load_final_clusters,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PrepareConfig:
    """Inputs, parameters, and parquet output for the fitness table preparation."""
    final_clusters: Path
    protein_features: Path
    grna_data: Path
    clip_upper: float
    output_fitness_table: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dir exists."""
        for path in [self.final_clusters, self.protein_features, self.grna_data]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_fitness_table.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PrepareConfig) -> None:
    """Load -> merge -> clip -> write the fitness table parquet intermediate."""
    config.validate()

    final_clusters = load_final_clusters(config.final_clusters)
    protein_features = pd.read_csv(config.protein_features, sep="\t")
    grna_data = pd.read_csv(config.grna_data, sep="\t")

    fitness_table = build_fitness_table(
        final_clusters, protein_features, grna_data, clip_upper=config.clip_upper
    )
    write_parquet(fitness_table, config.output_fitness_table)

    logger.success(f"Prepared fitness table: {len(fitness_table):,} genes")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare the merged fitness table parquet intermediate")
    parser.add_argument("--final-clusters", type=Path, required=True, help="Curated final_clusters.tsv")
    parser.add_argument("--protein-features", type=Path, required=True, help="pombe_coding_gene_protein_features.tsv")
    parser.add_argument("--grna-data", type=Path, required=True, help="gRNA HD-data fitted parameters TSV")
    parser.add_argument("--clip-upper", type=float, default=CLIP_UPPER, help="Upper cap for integration-density columns")
    parser.add_argument("--output-fitness-table", type=Path, required=True, help="Output fitness_table.parquet")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the preparation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PrepareConfig(
            final_clusters=args.final_clusters,
            protein_features=args.protein_features,
            grna_data=args.grna_data,
            clip_upper=args.clip_upper,
            output_fitness_table=args.output_fitness_table,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
