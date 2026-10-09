#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Comparison Figures
========================

Stage 3 of the comparison split: read the prepared fitness_table parquet
intermediate AND the fitness_correlation_stats.tsv, then render
- pairwise_fitness_comparison.pdf: ONE n x n scatter matrix (lower triangle,
  every pair of the fitness columns, variables ordered by similarity);
- correlation_heatmap.pdf: ONE figure holding the Pearson and Spearman
  correlation matrices side by side, with study-category colour bands.

Both figures take the order one average-linkage tree on 1 - r puts the columns
in (core.cluster_comparison_columns), read once and shared, so the two figures
always read the same way.

Output
------
- pairwise_fitness_comparison.pdf: n x n pairwise scatter matrix.
- correlation_heatmap.pdf: Pearson | Spearman correlation heatmaps.

Usage
-----
    python workflow/scripts/comparison/plot_comparison_figures.py \\
        --fitness-table results/6a_comparison/HD_DIT_HAP/_work/fitness_table.parquet \\
        --stats results/6a_comparison/HD_DIT_HAP/fitness_correlation_stats.tsv \\
        --output-dir results/6a_comparison/HD_DIT_HAP

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-10-08
Version:  4.0.0
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

matplotlib.use("Agg")  # headless: this script only writes figures, never displays
from loguru import logger  # noqa: E402

# 4. Local Imports (relative path resolution)
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402
from comparison.core import plot_comparison_figures, transform_fitness_columns  # noqa: E402


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True, slots=True)
class PlotConfig:
    """Parquet + stats TSV inputs and the output directory for the comparison figures."""
    fitness_table: Path
    stats: Path
    output_dir: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [self.fitness_table, self.stats]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_dir.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotConfig) -> None:
    """Read parquet + stats TSV -> plot the pairwise matrix + the correlation heatmap figure.

    The table is transformed the same way the stats rule transformed it, so the
    scatter matrix draws exactly the values the TSV's coefficients describe.
    """
    config.validate()

    fitness_table = transform_fitness_columns(read_parquet(config.fitness_table))
    stats = pd.read_csv(config.stats, sep="\t")
    columns = list(dict.fromkeys([*stats["col_x"], *stats["col_y"]]))

    plot_comparison_figures(fitness_table, stats, columns, config.output_dir)

    logger.success(f"Comparison figures: {len(stats):,} pairs plotted into {config.output_dir}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Pairwise scatter matrix + correlation heatmap figure")
    parser.add_argument("--fitness-table", type=Path, required=True, help="Input fitness_table.parquet")
    parser.add_argument("--stats", type=Path, required=True, help="Input fitness_correlation_stats.tsv")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory for figures")
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
