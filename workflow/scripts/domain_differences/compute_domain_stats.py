#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Intra-gene DR Heterogeneity (Domain-Difference) Candidate Statistics
====================================================================

Stage 2 of the domain-differences split: read the prepared gene_result /
annotations parquet intermediates, select genes with high gene-level DR
(> dr_threshold), restrict the annotations to in-gene insertions
(IN_GENE_FILTER), position each in-gene insertion along the CDS
(insertion_fraction), and aggregate into per-gene distribution statistics
(count, mean, std of insertion_fraction). Depends only on
prepare_domain_data's output.

See workflow/src/domain_differences/core.py for the full notebook-vs-script
deviation writeup (the source notebook is a visualization-only notebook that
computes no per-gene statistics table).

Output
------
- domain_candidate_stats.tsv: one row per high-DR gene that has >=1 in-gene
  insertion, columns [Systematic ID, Name, n_insertions,
  mean_insertion_fraction, std_insertion_fraction, gene_DR], sorted by
  std_insertion_fraction descending (single-insertion genes have NaN std and
  sort last).

Usage
-----
    python compute_domain_stats.py \\
        --gene-result .../_work/gene_result.parquet \\
        --annotations .../_work/annotations.parquet \\
        --dr-threshold 0.15 \\
        --output-stats results/domain_differences/{dataset}/domain_candidate_stats.tsv

Author:   Yusheng Yang (guidance) + Claude (implementation)
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

# 2. Third-party Imports
from loguru import logger

# 3. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io import read_parquet  # noqa: E402
from workflow.src.domain_differences.core import (  # noqa: E402
    DR_THRESHOLD,
    IN_GENE_FILTER,
    compute_domain_candidate_stats,
    compute_insertion_fraction,
    filter_high_dr_genes,
)


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class DomainConfig:
    """Parquet inputs, output, and threshold for the domain-difference candidate analysis."""
    gene_result: Path
    annotations: Path
    output_stats: Path
    dr_threshold: float = DR_THRESHOLD

    def validate(self) -> None:
        """Raise ValueError if any input is missing or the threshold is invalid, then ensure the output dir exists."""
        for path in [self.gene_result, self.annotations]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.dr_threshold < 0:
            raise ValueError(f"dr_threshold must be non-negative, got {self.dr_threshold}")
        self.output_stats.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: DomainConfig) -> None:
    """Read parquet -> select high-DR genes -> position in-gene insertions -> per-gene stats -> TSV."""
    config.validate()

    gene_result = read_parquet(config.gene_result)
    annotations = read_parquet(config.annotations)

    high_dr = filter_high_dr_genes(gene_result, config.dr_threshold)
    logger.info(f"Genes with gene-level DR > {config.dr_threshold}: {len(high_dr):,}")

    in_gene = annotations.query(IN_GENE_FILTER)
    logger.info(f"In-gene insertions (IN_GENE_FILTER): {len(in_gene):,}")
    in_gene = compute_insertion_fraction(in_gene)

    stats = compute_domain_candidate_stats(in_gene, high_dr)
    stats.to_csv(config.output_stats, sep="\t", index=False)

    n_multi = int((stats["n_insertions"] > 1).sum())
    logger.success(
        f"Domain candidates: {len(stats):,} high-DR genes with in-gene insertions "
        f"({n_multi:,} with >1 insertion), {int(stats['n_insertions'].sum()):,} insertions total "
        f"-> {config.output_stats}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Compute intra-gene DR heterogeneity (domain-difference) candidate statistics")
    parser.add_argument("--gene-result", type=Path, required=True, help="Prepared gene_result.parquet")
    parser.add_argument("--annotations", type=Path, required=True, help="Prepared annotations.parquet")
    parser.add_argument("--dr-threshold", type=float, default=DR_THRESHOLD, help="Gene-level DR selection cutoff (default: 0.15)")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output per-gene domain candidate stats TSV")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = DomainConfig(
            gene_result=args.gene_result,
            annotations=args.annotations,
            output_stats=args.output_stats,
            dr_threshold=args.dr_threshold,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
