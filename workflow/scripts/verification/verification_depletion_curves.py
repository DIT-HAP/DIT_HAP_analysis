#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Verification Depletion Curves
=============================

Stage 2c of the verification split: read the prepared merged parquet plus the
gene-level per-timepoint statistics (and, for HD, the curated gRNA fitted
parameters), and emit one single-page depletion-curve figure per critical-gene
group. Each panel is a gene's DIT-HAP measured points + Gompertz fit + inflection
tangent, with the gRNA curve overlaid where the dataset has one. Genes per group
are the same filter-selected outliers the boxplot rule uses. Depends only on
prepare_verification_table's merged.parquet plus the timepoint inputs.

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
from figure_render.verification import render_depletion_curves_figure  # noqa: E402
from logging_setup import setup_logger  # noqa: E402
from release_schema import read_gene_level  # noqa: E402
from verification.core import (  # noqa: E402
    CRITICAL_GROUPS,
    load_grna_timepoints,
    select_group_outliers,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class DepletionCurveConfig:
    """Merged parquet + timepoint inputs, figure output dir. grna_timepoints is optional (HD-only)."""
    merged: Path
    gene_timepoints: Path
    output_dir: Path
    grna_timepoints: Path | None = None

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [self.merged, self.gene_timepoints]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.grna_timepoints is not None and not self.grna_timepoints.exists():
            raise ValueError(f"gRNA timepoints given but not found: {self.grna_timepoints}")
        self.output_dir.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: DepletionCurveConfig) -> None:
    """Read merged parquet + timepoints -> write one depletion-curve figure per critical group."""
    config.validate()

    merged = read_parquet(config.merged)
    gene_timepoints = read_gene_level(config.gene_timepoints).set_index("Systematic ID")
    grna_timepoints = load_grna_timepoints(config.grna_timepoints)

    for group in CRITICAL_GROUPS:
        render_depletion_curves_figure(
            select_group_outliers(merged, group),
            gene_timepoints,
            grna_timepoints,
            config.output_dir / f"depletion_curves_{group}",
            group=group,
        )

    logger.success(f"Depletion curves written to {config.output_dir}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Verification DIT-HAP (+gRNA) depletion curves per critical group")
    parser.add_argument("--merged", type=Path, required=True, help="Input merged.parquet")
    parser.add_argument("--gene-timepoints", type=Path, required=True, help="Gene-level fitting statistics with YES0-4")
    parser.add_argument("--grna-timepoints", type=Path, default=None, help="Curated gRNA fitted-parameters TSV (optional; HD-only overlay)")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output dir for per-group depletion-curve figures")
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
            output_dir=args.output_dir,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
