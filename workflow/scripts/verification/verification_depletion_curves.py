#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Verification Depletion Curves
=============================

Stage 2c of the verification split: read the prepared merged parquet + the
gene-level per-timepoint statistics (and, for HD, the gRNA per-timepoint LFC),
and emit the depletion-curve PDF — one 4-column grid per critical-gene group,
each panel a gene's DIT-HAP measured points + Gompertz fit + inflection slope,
with the gRNA curve overlaid where available. Genes per group are the same
filter-selected outliers the boxplot rule uses. Depends only on
prepare_verification_table's merged.parquet plus the timepoint inputs.

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
    build_depletion_curve_pdf,
    load_grna_timepoints,
    load_gene_level,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class DepletionCurveConfig:
    """Merged parquet + timepoint inputs, PDF output. grna_timepoints is optional (HD-only)."""
    merged: Path
    gene_timepoints: Path
    output_depletion_curves: Path
    grna_timepoints: Path | None = None

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dir exists."""
        for path in [self.merged, self.gene_timepoints]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.grna_timepoints is not None and not self.grna_timepoints.exists():
            raise ValueError(f"gRNA timepoints given but not found: {self.grna_timepoints}")
        self.output_depletion_curves.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: DepletionCurveConfig) -> None:
    """Read merged parquet + timepoints -> write the depletion-curve PDF."""
    config.validate()

    merged = read_parquet(config.merged)
    gene_timepoints = load_gene_level(config.gene_timepoints).set_index("Systematic ID")
    grna_timepoints = load_grna_timepoints(config.grna_timepoints)

    build_depletion_curve_pdf(merged, gene_timepoints, grna_timepoints, config.output_depletion_curves)

    logger.success(f"Depletion curves written to {config.output_depletion_curves}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Verification DIT-HAP (+gRNA) depletion curves per critical group")
    parser.add_argument("--merged", type=Path, required=True, help="Input merged.parquet")
    parser.add_argument("--gene-timepoints", type=Path, required=True, help="Gene-level fitting statistics with YES0-4")
    parser.add_argument("--grna-timepoints", type=Path, default=None, help="gRNA per-timepoint LFC CSV (optional; HD-only overlay)")
    parser.add_argument("--output-depletion-curves", type=Path, required=True, help="Output depletion-curve PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the depletion-curve stage, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = DepletionCurveConfig(
            merged=args.merged,
            gene_timepoints=args.gene_timepoints,
            grna_timepoints=args.grna_timepoints,
            output_depletion_curves=args.output_depletion_curves,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
