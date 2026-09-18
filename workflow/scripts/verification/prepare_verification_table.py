#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare Verification Tables
===========================

Stage 1 of the verification split (see
docs/plans/2026-07-22-verification-rules-split-design.md): load the gene-level
DIT-HAP results, deletion-library categories, and curated essentiality
verification table, then merge them into the two parquet intermediates consumed
by the category-summary / boxplot / depletion-curve rules:

- merged.parquet: gene-level DR/DL + DeletionLibrary_essentiality + Category
  (raw curated label) + Category_with_essentiality (one row per gene).
- verification.parquet: the curated verification table, all columns, plus the
  simplified `Verification result` the outlier bucketing groups on.

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
from io_table import write_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402
from release_schema import read_gene_level  # noqa: E402
from verification.core import (  # noqa: E402
    CATEGORY_WITH_ESSENTIALITY_COLUMN,
    apply_category_with_essentiality,
    load_deletion_library,
    load_verification,
    merge_deletion_library,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PrepareConfig:
    """Inputs and parquet outputs for the verification table preparation."""
    fitting_results: Path
    deletion_library: Path
    essentiality_verification: Path
    output_merged: Path
    output_verification: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.fitting_results, self.deletion_library, self.essentiality_verification]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_merged, self.output_verification]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PrepareConfig) -> None:
    """Load -> merge -> write the two parquet intermediates."""
    config.validate()

    merged = merge_deletion_library(
        read_gene_level(config.fitting_results),
        load_deletion_library(config.deletion_library),
    )
    merged[CATEGORY_WITH_ESSENTIALITY_COLUMN] = merged.apply(apply_category_with_essentiality, axis=1)
    verification = load_verification(config.essentiality_verification)

    write_parquet(merged, config.output_merged)
    write_parquet(verification, config.output_verification)

    logger.success(
        f"Prepared verification tables: {len(merged):,} genes, "
        f"{len(verification):,} verified-gene rows"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare verification parquet intermediates")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Gene-level fitting_results.tsv")
    parser.add_argument("--deletion-library", type=Path, required=True, help="Curated deletion_library_categories.xlsx")
    parser.add_argument("--essentiality-verification", type=Path, required=True, help="Curated essentiality_verification.csv")
    parser.add_argument("--output-merged", type=Path, required=True, help="Output merged.parquet")
    parser.add_argument("--output-verification", type=Path, required=True, help="Output verification.parquet")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the preparation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PrepareConfig(
            fitting_results=args.fitting_results,
            deletion_library=args.deletion_library,
            essentiality_verification=args.essentiality_verification,
            output_merged=args.output_merged,
            output_verification=args.output_verification,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
