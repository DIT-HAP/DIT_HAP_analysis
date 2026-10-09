#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
QC Fitness Distributions
========================

Distribution QC for the comparison stage: every fitness column the comparison
correlates is drawn twice -- raw and log10-transformed -- so the choice of
statistic can be checked against the data.

A column whose two panels do not look alike is heavy-tailed, and its Pearson r
is carried by the tail where its Spearman rho is not. Non-positive values are
dropped from the log10 panels (and counted in the log) rather than shifted.

Rows run in the order both comparison figures use.

Output
------
- fitness_distributions_qc.pdf: one row of panels per column (raw | log10).
- fitness_distributions_qc.review.png: the same figure, for review.

Usage
-----
    python workflow/scripts/comparison/qc_fitness_distributions.py \\
        --fitness-table results/6a_comparison/HD_DIT_HAP/_work/fitness_table.parquet \\
        --stats results/6a_comparison/HD_DIT_HAP/fitness_correlation_stats.tsv \\
        --output-dir results/6a_comparison/HD_DIT_HAP

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-10-09
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

matplotlib.use("Agg")  # headless: this script only writes figures, never displays
from loguru import logger  # noqa: E402

# 4. Local Imports (relative path resolution)
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402
from comparison.core import (  # noqa: E402
    cluster_comparison_columns,
    plot_fitness_distributions,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True, slots=True)
class QcConfig:
    """Parquet + stats TSV inputs and the output directory for the distribution QC."""
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
def run(config: QcConfig) -> None:
    """Read parquet + stats TSV -> plot every correlated column's raw and log10 distribution."""
    config.validate()

    fitness_table = read_parquet(config.fitness_table)
    stats = pd.read_csv(config.stats, sep="\t")
    columns = list(dict.fromkeys([*stats["col_x"], *stats["col_y"]]))

    plot_fitness_distributions(
        fitness_table,
        config.output_dir / "fitness_distributions_qc",
        order=cluster_comparison_columns(stats, columns)[0],
    )

    logger.success(f"Fitness distribution QC: {len(columns)} columns into {config.output_dir}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Raw and log10 distribution QC per fitness column")
    parser.add_argument("--fitness-table", type=Path, required=True, help="Input fitness_table.parquet")
    parser.add_argument("--stats", type=Path, required=True, help="Input fitness_correlation_stats.tsv")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory for the QC figure")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the QC plotting, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = QcConfig(
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
