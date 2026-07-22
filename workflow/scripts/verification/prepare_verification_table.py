#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare Verification Tables
===========================

Stage 1 of the verification split (see
docs/plans/2026-07-22-verification-rules-split-design.md): load the gene-level
DIT-HAP results, deletion-library categories, and curated essentiality
verification table, then merge them into three parquet intermediates consumed
by the category-summary / boxplot / depletion-curve rules:

- merged.parquet: gene-level DR/DL + DeletionLibrary_essentiality + Category +
  Category_with_essentiality + cat_canon (one row per gene).
- final_merged.parquet: the curated-verification genes with area day3-6 columns
  (feeds the critical-gene review TSVs).
- simplified_verification.parquet: Systematic ID / Verification result /
  Verified essentiality (+ the manual gpd1 row), for outlier bucketing.

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
from workflow.src.io import write_parquet  # noqa: E402
from workflow.src.verification.core import (  # noqa: E402
    apply_category_with_essentiality,
    build_final_merged,
    canonicalize_category,
    load_deletion_library,
    load_essentiality_verification,
    load_essentiality_verification_full,
    load_gene_level,
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
    output_final_merged: Path
    output_simplified_verification: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.fitting_results, self.deletion_library, self.essentiality_verification]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_merged, self.output_final_merged, self.output_simplified_verification]:
            out.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PrepareConfig) -> None:
    """Load -> merge -> write the three parquet intermediates."""
    config.validate()

    gene_result = load_gene_level(config.fitting_results)
    deletion_library = load_deletion_library(config.deletion_library)
    simplified_verification = load_essentiality_verification(config.essentiality_verification)
    verification_full = load_essentiality_verification_full(config.essentiality_verification)

    merged = merge_deletion_library(gene_result, deletion_library)
    merged["Category_with_essentiality"] = merged.apply(apply_category_with_essentiality, axis=1)
    merged = canonicalize_category(merged)
    final_merged = build_final_merged(merged, verification_full)

    write_parquet(merged, config.output_merged)
    write_parquet(final_merged, config.output_final_merged)
    write_parquet(simplified_verification, config.output_simplified_verification)

    logger.success(
        f"Prepared verification tables: {len(merged):,} genes, "
        f"{len(final_merged):,} verified-gene rows, {len(simplified_verification):,} simplified verifications"
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
    parser.add_argument("--output-final-merged", type=Path, required=True, help="Output final_merged.parquet")
    parser.add_argument("--output-simplified-verification", type=Path, required=True, help="Output simplified_verification.parquet")
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
            output_final_merged=args.output_final_merged,
            output_simplified_verification=args.output_simplified_verification,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
