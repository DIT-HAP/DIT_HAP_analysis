#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute Coverage Stats
========================

Stage 2a of the coverage split: read the prepared annotations /
gene_result parquet intermediates and emit the coverage_stats.tsv (insertion
+ gene + essential + non_essential + per-chromosome coverage). Depends only
on prepare_coverage_data's output, so it re-runs independently of the figures
rule.

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

# 2. Third-party Imports
from loguru import logger

# 3. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io import read_parquet  # noqa: E402
from workflow.src.coverage.core import (  # noqa: E402
    build_detailed_gene_table,
    build_stats_table,
    compute_characterisation_status_coverage,
    compute_deletion_viability_coverage,
    compute_essentiality_category_coverage,
    compute_essentiality_coverage,
    compute_gene_coverage,
    compute_insertion_coverage,
    compute_per_chromosome_insertion_coverage,
    write_detailed_gene_excel,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class ComputeStatsConfig:
    """Parquet inputs + TSV output for the coverage stats computation."""
    annotations: Path
    gene_result: Path
    gene_metadata: Path
    output_stats: Path
    output_detailed_genes_xlsx: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.annotations, self.gene_result, self.gene_metadata]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_stats.parent.mkdir(parents=True, exist_ok=True)
        self.output_detailed_genes_xlsx.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: ComputeStatsConfig) -> None:
    """Read parquet -> compute coverage stats -> write TSV."""
    config.validate()

    annotations = read_parquet(config.annotations)
    gene_result = read_parquet(config.gene_result)
    gene_metadata = read_parquet(config.gene_metadata)

    insertion_coverage = compute_insertion_coverage(annotations)
    gene_coverage = compute_gene_coverage(gene_result)
    essentiality_coverage = compute_essentiality_coverage(gene_result)
    per_chromosome = compute_per_chromosome_insertion_coverage(annotations)
    characterisation_status_coverage = compute_characterisation_status_coverage(gene_result)
    deletion_viability_coverage = compute_deletion_viability_coverage(gene_result)
    essentiality_category_coverage = compute_essentiality_category_coverage(gene_result)

    stats_table = build_stats_table(
        insertion_coverage,
        gene_coverage,
        essentiality_coverage,
        per_chromosome,
        characterisation_status_coverage,
        deletion_viability_coverage,
        essentiality_category_coverage,
    )
    stats_table.to_csv(config.output_stats, sep="\t", index=False)

    # Build and write detailed gene-level Excel with multiple sheets
    detailed_genes = build_detailed_gene_table(gene_result, gene_metadata)
    write_detailed_gene_excel(detailed_genes, config.output_detailed_genes_xlsx)

    logger.success(
        f"Coverage: {insertion_coverage['in_gene']:,}/{insertion_coverage['total']:,} insertions in-gene, "
        f"{gene_coverage['covered']:,}/{gene_coverage['total']:,} genes covered "
        f"({essentiality_coverage['essential']['covered']:,}/{essentiality_coverage['essential']['total']:,} essential)"
    )
    if characterisation_status_coverage:
        logger.info(f"Computed coverage for {len(characterisation_status_coverage)} characterisation_status categories")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Compute gene insertion coverage statistics")
    parser.add_argument("--annotations", type=Path, required=True, help="Input annotations.parquet")
    parser.add_argument("--gene-result", type=Path, required=True, help="Input gene_result.parquet")
    parser.add_argument("--gene-metadata", type=Path, required=True, help="Input gene metadata parquet")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output coverage stats TSV")
    parser.add_argument("--output-detailed-genes-xlsx", type=Path, required=True, help="Output detailed gene Excel file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the computation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = ComputeStatsConfig(
            annotations=args.annotations,
            gene_result=args.gene_result,
            gene_metadata=args.gene_metadata,
            output_stats=args.output_stats,
            output_detailed_genes_xlsx=args.output_detailed_genes_xlsx,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
