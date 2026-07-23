#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Coverage Figures
=======================

Stage 2b of the coverage split: read the prepared annotations / gene_result
parquet intermediates, recompute the small coverage dicts (same core
functions as compute_coverage_stats.py, for symmetry), and emit the
coverage_figures.pdf (donut charts + per-chromosome bars + DR/DL histograms).
Depends only on prepare_coverage_data's output, so it re-runs independently
of the stats-table rule.

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
import matplotlib

matplotlib.use("Agg")  # headless: this script only writes a PDF, never displays
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from loguru import logger  # noqa: E402

# 3. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io import read_parquet  # noqa: E402
from workflow.src.coverage.core import (  # noqa: E402
    compute_essentiality_coverage,
    compute_gene_coverage,
    compute_insertion_coverage,
    compute_per_chromosome_insertion_coverage,
    plot_coverage_donuts,
    plot_dr_dl_histograms,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PlotFiguresConfig:
    """Parquet inputs + PDF output for the coverage figures."""
    annotations: Path
    gene_result: Path
    output_figures: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.annotations, self.gene_result]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_figures.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotFiguresConfig) -> None:
    """Read parquet -> recompute coverage dicts -> write figures PDF."""
    config.validate()

    annotations = read_parquet(config.annotations)
    gene_result = read_parquet(config.gene_result)

    insertion_coverage = compute_insertion_coverage(annotations)
    gene_coverage = compute_gene_coverage(gene_result)
    essentiality_coverage = compute_essentiality_coverage(gene_result)
    per_chromosome = compute_per_chromosome_insertion_coverage(annotations)

    fig_donuts = plot_coverage_donuts(insertion_coverage, gene_coverage, essentiality_coverage, per_chromosome)
    fig_hist = plot_dr_dl_histograms(gene_result)

    with PdfPages(config.output_figures) as pdf:
        pdf.savefig(fig_donuts, dpi=300, bbox_inches="tight")
        pdf.savefig(fig_hist, dpi=300, bbox_inches="tight")
    plt.close(fig_donuts)
    plt.close(fig_hist)

    logger.success(f"Wrote coverage figures: {config.output_figures}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Plot gene insertion coverage figures")
    parser.add_argument("--annotations", type=Path, required=True, help="Input annotations.parquet")
    parser.add_argument("--gene-result", type=Path, required=True, help="Input gene_result.parquet")
    parser.add_argument("--output-figures", type=Path, required=True, help="Output coverage figures PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the plotting, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotFiguresConfig(
            annotations=args.annotations,
            gene_result=args.gene_result,
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
