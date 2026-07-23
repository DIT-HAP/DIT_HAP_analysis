#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Comparison Figures
=========================

Stage 3 of the comparison split: read the prepared fitness_table parquet
intermediate AND the fitness_correlation_stats.tsv (for the col_x/col_y pairs
that SURVIVED the per-pair overlap filter in compute_comparison_stats), then
render the pairwise scatter-matrix PDF with a Gaussian-KDE density overlay per
panel. Driving the plot grid from the stats TSV's surviving pairs (rather than
recomputing them) keeps the PDF panels and TSV rows in permanent agreement
even though the two rules now run independently.

Output
------
- pairwise_fitness_comparison.pdf: scatter matrix (one panel per pair) with a
  Gaussian-KDE density overlay and Pearson r/p/n annotation.

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

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
import matplotlib

matplotlib.use("Agg")  # headless: this script only writes a PDF, never displays
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io import read_parquet  # noqa: E402
from workflow.src.comparison.core import plot_pairwise_comparison  # noqa: E402


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PlotConfig:
    """Parquet + stats TSV inputs and PDF output for the pairwise comparison figures."""
    fitness_table: Path
    stats: Path
    output_figures: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [self.fitness_table, self.stats]:
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
def run(config: PlotConfig) -> None:
    """Read parquet + stats TSV -> plot the surviving pairs -> write the PDF."""
    config.validate()

    fitness_table = read_parquet(config.fitness_table)
    stats = pd.read_csv(config.stats, sep="\t")

    # Drive the plot grid from the surviving stats pairs (post per-pair overlap
    # filter) so PDF panels and TSV rows always agree on which pairs exist.
    surviving_pairs = list(zip(stats["col_x"], stats["col_y"]))
    fig = plot_pairwise_comparison(fitness_table, surviving_pairs)
    with PdfPages(config.output_figures) as pdf:
        pdf.savefig(fig, dpi=300, bbox_inches="tight")
    plt.close(fig)

    logger.success(f"Comparison figures: {len(surviving_pairs):,} panels plotted")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Pairwise fitness comparison scatter-matrix PDF")
    parser.add_argument("--fitness-table", type=Path, required=True, help="Input fitness_table.parquet")
    parser.add_argument("--stats", type=Path, required=True, help="Input fitness_correlation_stats.tsv")
    parser.add_argument("--output-figures", type=Path, required=True, help="Output pairwise comparison PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the plotting, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotConfig(
            fitness_table=args.fitness_table,
            stats=args.stats,
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
