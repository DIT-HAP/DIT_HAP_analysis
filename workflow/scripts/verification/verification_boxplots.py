#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Verification Boxplots
=====================

Stage 2b of the verification split: read the prepared merged / verification
parquet intermediates and emit one figure per analysis — the canonical-category
DR violin plot, plus a DR-violin-versus-verification-donut figure for each
critical-gene group — and the per-group critical-gene review TSVs. Depends only
on prepare_verification_table's output.

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
from figure_render.verification import (  # noqa: E402
    render_category_boxplot_figure,
    render_critical_group_figure,
)
from logging_setup import setup_logger  # noqa: E402
from verification.core import (  # noqa: E402
    CRITICAL_GROUPS,
    build_final_merged,
    critical_group_boxplot_data,
    order_verification_buckets,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class BoxplotConfig:
    """Parquet inputs + figure/TSV output paths for the boxplot stage."""
    merged: Path
    verification: Path
    output_figure: Path
    output_critical_genes_dir: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.merged, self.verification]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_figure.parent.mkdir(parents=True, exist_ok=True)
        self.output_critical_genes_dir.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: BoxplotConfig) -> None:
    """Read parquet -> write the category box plot, per-group figures and review TSVs."""
    config.validate()

    merged = read_parquet(config.merged)
    verification = read_parquet(config.verification)
    # The review TSVs carry each gene's metrics and colony-area columns together,
    # so the join happens here rather than being written out as a third
    # intermediate that nothing else reads.
    final_merged = build_final_merged(merged, verification)

    render_category_boxplot_figure(merged, config.output_figure)

    order = order_verification_buckets(verification["Verification result"].dropna())
    for group in CRITICAL_GROUPS:
        dr_by_bucket, detail = critical_group_boxplot_data(
            merged, final_merged, verification, group
        )
        detail.to_csv(
            config.output_critical_genes_dir / f"critical_genes_{group}.tsv",
            sep="\t",
            index=False,
        )
        render_critical_group_figure(
            dr_by_bucket,
            config.output_critical_genes_dir / f"critical_genes_{group}",
            group=group,
            order=order,
        )
        logger.info(f"{group}: {len(detail):,} genes across {len(dr_by_bucket)} verification buckets")

    logger.success(f"Box plot figures + critical-gene TSVs written to {config.output_critical_genes_dir.parent}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Verification category box plot + per-critical-group figures and review TSVs")
    parser.add_argument("--merged", type=Path, required=True, help="Input merged.parquet")
    parser.add_argument("--verification", type=Path, required=True, help="Input verification.parquet")
    parser.add_argument("--output-figure", type=Path, required=True, help="Output stem for the category box plot (writes <stem>.pdf + <stem>.review.png)")
    parser.add_argument("--output-critical-genes-dir", type=Path, required=True, help="Output dir for per-group figures + review TSVs")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the boxplot stage, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = BoxplotConfig(
            merged=args.merged,
            verification=args.verification,
            output_figure=args.output_figure,
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
