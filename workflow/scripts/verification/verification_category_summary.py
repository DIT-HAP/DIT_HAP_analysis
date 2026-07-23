#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Verification Category Summary
=============================

Stage 2a of the verification split: read the prepared merged /
simplified_verification parquet intermediates and emit the category-level
stats TSV + the two-page summary PDF (phenotype-category donut + DR-by-category
scatter). Depends only on prepare_verification_table's output, so it re-runs
independently of the boxplot / depletion-curve rules.

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
from workflow.src.verification.core import (  # noqa: E402
    build_category_summary_pdf,
    build_stats_table,
    compute_category_stats,
    compute_category_with_essentiality_stats,
    compute_verification_match_stats,
    merge_essentiality_verification,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class CategorySummaryConfig:
    """Parquet inputs + TSV/PDF outputs for the category summary."""
    merged: Path
    simplified_verification: Path
    output_stats: Path
    output_figures: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.merged, self.simplified_verification]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_stats, self.output_figures]:
            out.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: CategorySummaryConfig) -> None:
    """Read parquet -> compute stats -> write TSV + summary PDF."""
    config.validate()

    merged = read_parquet(config.merged)
    simplified_verification = read_parquet(config.simplified_verification)
    merged_with_verification = merge_essentiality_verification(merged, simplified_verification)

    category_stats = compute_category_stats(merged)
    category_with_essentiality_stats = compute_category_with_essentiality_stats(merged)
    verification_stats = compute_verification_match_stats(merged_with_verification)

    stats_table = build_stats_table(category_stats, category_with_essentiality_stats, verification_stats)
    stats_table.to_csv(config.output_stats, sep="\t", index=False)

    build_category_summary_pdf(category_stats, merged, config.output_figures)

    logger.success(
        f"Category summary: {len(merged):,} genes across {len(category_stats):,} categories, "
        f"{verification_stats['match']:,}/{verification_stats['verified_total']:,} curated verifications match"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Verification category summary (stats TSV + donut/scatter PDF)")
    parser.add_argument("--merged", type=Path, required=True, help="Input merged.parquet")
    parser.add_argument("--simplified-verification", type=Path, required=True, help="Input simplified_verification.parquet")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output verification stats TSV")
    parser.add_argument("--output-figures", type=Path, required=True, help="Output category summary PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the summary, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = CategorySummaryConfig(
            merged=args.merged,
            simplified_verification=args.simplified_verification,
            output_stats=args.output_stats,
            output_figures=args.output_figures,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
