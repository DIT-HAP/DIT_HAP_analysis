#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Verification Category Summary
=============================

Stage 2a of the verification split: read the prepared merged / verification
parquet intermediates and emit the category-level stats TSV plus the
deletion-library comparison figure (phenotype-category donut + DR-by-category
scatter). Depends only on prepare_verification_table's output,
so it re-runs independently of the boxplot / depletion-curve rules.

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-22
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
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402
from figure_render.verification import render_category_summary_figure  # noqa: E402
from verification.core import (  # noqa: E402
    CATEGORY_COLUMN,
    build_stats_table,
    compute_verification_match_stats,
    order_categories,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class CategorySummaryConfig:
    """Parquet inputs + TSV/summary-figure outputs for the category summary."""
    merged: Path
    verification: Path
    output_stats: Path
    output_figure: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.merged, self.verification]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_stats, self.output_figure]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: CategorySummaryConfig) -> None:
    """Read parquet -> compute stats -> write TSV + comparison figure."""
    config.validate()

    merged = read_parquet(config.merged)
    verification = read_parquet(config.verification)

    build_stats_table(merged, verification).to_csv(config.output_stats, sep="\t", index=False)

    render_category_summary_figure(
        merged,
        config.output_figure,
        order=order_categories(merged[CATEGORY_COLUMN]),
    )

    match_stats = compute_verification_match_stats(merged, verification)
    logger.success(
        f"Category summary: {len(merged):,} genes across "
        f"{merged[CATEGORY_COLUMN].nunique():,} categories, "
        f"{match_stats['match']:,}/{match_stats['verified_total']:,} curated verifications match"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Verification category summary (stats TSV + comparison figure)")
    parser.add_argument("--merged", type=Path, required=True, help="Input merged.parquet")
    parser.add_argument("--verification", type=Path, required=True, help="Input verification.parquet")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output verification stats TSV")
    parser.add_argument("--output-figure", type=Path, required=True, help="Output figure stem (writes <stem>.pdf + <stem>.review.png)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the summary, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = CategorySummaryConfig(
            merged=args.merged,
            verification=args.verification,
            output_stats=args.output_stats,
            output_figure=args.output_figure,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
