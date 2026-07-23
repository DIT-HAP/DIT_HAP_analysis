#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute Comparison Stats
=========================

Stage 2 of the comparison split: read the prepared fitness_table parquet
intermediate, select the fitness columns with enough non-NaN data, and emit
the long-form Pearson r/p_value/n stats TSV. Depends only on
prepare_fitness_table's output, so it re-runs independently of the figures
rule (which reads this stats TSV back to know which pairs survived).

Output
------
- fitness_correlation_stats.tsv: long-form (col_x, col_y, pair, r, p_value, n),
  one row per unordered pair of available fitness columns.

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
from workflow.src.comparison.core import (  # noqa: E402
    compute_correlation_stats,
    select_fitness_columns,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class StatsConfig:
    """Parquet input + TSV output for the correlation stats computation."""
    fitness_table: Path
    output_stats: Path

    def validate(self) -> None:
        """Raise ValueError if the required input is missing, then ensure the output dir exists."""
        if not self.fitness_table.exists():
            raise ValueError(f"Required input not found: {self.fitness_table}")
        self.output_stats.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: StatsConfig) -> None:
    """Read parquet -> select columns -> compute correlation stats -> write TSV."""
    config.validate()

    fitness_table = read_parquet(config.fitness_table)
    columns = select_fitness_columns(fitness_table)

    stats = compute_correlation_stats(fitness_table, columns)
    stats.to_csv(config.output_stats, sep="\t", index=False)

    logger.success(
        f"Comparison stats: {len(stats):,} fitness-column pairs correlated across "
        f"{len(columns):,} columns ({len(fitness_table):,} genes)"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Pairwise fitness correlation stats TSV")
    parser.add_argument("--fitness-table", type=Path, required=True, help="Input fitness_table.parquet")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output correlation stats TSV")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the stats computation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = StatsConfig(
            fitness_table=args.fitness_table,
            output_stats=args.output_stats,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
