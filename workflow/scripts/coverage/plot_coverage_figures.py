#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Coverage Figures
=======================

Stage 2b of the coverage split: emit coverage_figures.pdf. Donut charts +
per-chromosome bars read straight from coverage_stats.tsv (the numbers
compute_coverage_stats wrote), so the figures can never disagree with the
table. The DR/DL histograms still read the gene_result parquet, because they
need the per-gene DR/DL values that the aggregated stats table doesn't carry.

Pages:
  1. Overall coverage donuts (insertion/gene/essential/non-essential) + per-chromosome bars
  2. Overall DR/DL histograms (all/essential/non-essential rows)
  3. Per-characterisation_status coverage donuts (one per category)
  4. Per-characterisation_status DR/DL histograms (one row per category)

Depending on both the stats TSV and the gene_result parquet means editing the
stats rule now does force the figures to rebuild — the deliberate trade for
guaranteed figure/table agreement.

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
import matplotlib

matplotlib.use("Agg")  # headless: this script only writes a PDF, never displays
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from loguru import logger  # noqa: E402

# 3. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import pandas as pd  # noqa: E402
from workflow.src.io import read_parquet  # noqa: E402
from workflow.src.coverage.core import (  # noqa: E402
    coverage_dicts_from_stats_table,
    plot_characterisation_status_donuts,
    plot_characterisation_status_histograms,
    plot_coverage_donuts,
    plot_dr_dl_histograms,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PlotFiguresConfig:
    """Inputs (stats TSV + gene_result parquet) + PDF output for the coverage figures."""
    stats: Path
    gene_result: Path
    output_figures: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.stats, self.gene_result]:
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
    """Read stats TSV (donuts) + gene_result parquet (histograms) -> write figures PDF."""
    config.validate()

    stats = pd.read_csv(config.stats, sep="\t")
    gene_result = read_parquet(config.gene_result)

    # Donuts + per-chromosome bars come straight from the stats table, so the
    # figures report exactly what compute_coverage_stats wrote (no re-derivation).
    (
        insertion_coverage,
        gene_coverage,
        essentiality_coverage,
        per_chromosome,
        characterisation_status_coverage,
    ) = coverage_dicts_from_stats_table(stats)

    fig_donuts = plot_coverage_donuts(insertion_coverage, gene_coverage, essentiality_coverage, per_chromosome)
    fig_hist = plot_dr_dl_histograms(gene_result)

    figures = [fig_donuts, fig_hist]

    # Per-characterisation_status pages: donuts (from stats) + DR/DL histograms (from parquet).
    if characterisation_status_coverage:
        figures.append(plot_characterisation_status_donuts(characterisation_status_coverage))
        figures.append(plot_characterisation_status_histograms(gene_result))
        logger.info(f"Added per-category figures for {len(characterisation_status_coverage)} characterisation_status categories")

    with PdfPages(config.output_figures) as pdf:
        for fig in figures:
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
    for fig in figures:
        plt.close(fig)

    logger.success(f"Wrote coverage figures ({len(figures)} pages): {config.output_figures}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Plot gene insertion coverage figures")
    parser.add_argument("--stats", type=Path, required=True, help="Input coverage_stats.tsv (donuts read from here)")
    parser.add_argument("--gene-result", type=Path, required=True, help="Input gene_result.parquet (histograms read from here)")
    parser.add_argument("--output-figures", type=Path, required=True, help="Output coverage figures PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the plotting, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotFiguresConfig(
            stats=args.stats,
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
