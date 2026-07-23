#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Verification Boxplots
=====================

Stage 2b of the verification split: read the prepared merged / final_merged /
simplified_verification parquet intermediates and emit the boxplot+violin PDF
(basic per-category DR + four critical-gene groups, each with a boxplot and a
verification-composition donut) plus the per-group critical-gene review TSVs.
Depends only on prepare_verification_table's output.

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
from workflow.src.verification.core import build_boxplot_pdf  # noqa: E402


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class BoxplotConfig:
    """Parquet inputs + PDF/TSV outputs for the boxplot stage."""
    merged: Path
    final_merged: Path
    simplified_verification: Path
    output_boxplots: Path
    output_critical_genes_dir: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.merged, self.final_merged, self.simplified_verification]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_boxplots.parent.mkdir(parents=True, exist_ok=True)
        self.output_critical_genes_dir.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: BoxplotConfig) -> None:
    """Read parquet -> write boxplot/donut PDF + critical-gene review TSVs."""
    config.validate()

    merged = read_parquet(config.merged)
    final_merged = read_parquet(config.final_merged)
    simplified_verification = read_parquet(config.simplified_verification)

    build_boxplot_pdf(
        merged, final_merged, simplified_verification,
        config.output_boxplots, config.output_critical_genes_dir,
    )

    logger.success(f"Boxplots + critical-gene TSVs written to {config.output_boxplots.parent}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Verification boxplot/violin + critical-gene donuts + review TSVs")
    parser.add_argument("--merged", type=Path, required=True, help="Input merged.parquet")
    parser.add_argument("--final-merged", type=Path, required=True, help="Input final_merged.parquet")
    parser.add_argument("--simplified-verification", type=Path, required=True, help="Input simplified_verification.parquet")
    parser.add_argument("--output-boxplots", type=Path, required=True, help="Output boxplot/violin + donut PDF")
    parser.add_argument("--output-critical-genes-dir", type=Path, required=True, help="Output dir for per-group critical-gene TSVs")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the boxplot stage, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = BoxplotConfig(
            merged=args.merged,
            final_merged=args.final_merged,
            simplified_verification=args.simplified_verification,
            output_boxplots=args.output_boxplots,
            output_critical_genes_dir=args.output_critical_genes_dir,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
