#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Spike-In Correlation
============================

Stage 2b of the spike-in split: read the prepared spike_in_stats parquet
intermediate, drop the non-finite Relative_Read_Ratio row (the lowest-read
sample per site floors at Reads=0 -> log2(0) = -inf), fit a log-log linear
regression of read ratio vs dilution ratio, and emit the correlation PDF.
Depends only on prepare_spikein_data's output, so it re-runs independently of
the stats-TSV rule.

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
import numpy as np

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io import read_parquet  # noqa: E402
from workflow.src.spikein.core import compute_linear_regression_stats, plot_spike_in_correlation  # noqa: E402


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PlotCorrelationConfig:
    """Parquet input + PDF output for the spike-in correlation figure."""
    spike_in_stats: Path
    output_figure: Path

    def validate(self) -> None:
        """Raise ValueError if the input is missing, then ensure the output dir exists."""
        if not self.spike_in_stats.exists():
            raise ValueError(f"Required input not found: {self.spike_in_stats}")
        self.output_figure.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotCorrelationConfig) -> None:
    """Read parquet -> filter to finite rows -> fit -> save the correlation PDF."""
    config.validate()

    spike_in_stats = read_parquet(config.spike_in_stats)

    # The lowest-read sample per site floors at Reads=0 -> Relative_Read_Ratio
    # = log2(0) = -inf (kept in the TSV as the true computed value). Excluded
    # here so the fit/plot aren't skewed or broken by a non-finite point.
    finite = spike_in_stats[np.isfinite(spike_in_stats["Relative_Read_Ratio"])]

    stats = compute_linear_regression_stats(
        finite["Relative_Dilution_Ratio"], finite["Relative_Read_Ratio"]
    )
    plot_spike_in_correlation(finite, stats, config.output_figure)

    logger.success(
        f"Spike-in linearity: slope={stats['slope']:.3f}, R2={stats['r2']:.3f}, "
        f"p={stats['p_value']:.2e} ({len(finite)}/{len(spike_in_stats)} finite site x sample rows)"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Spike-in correlation figure (from the prepared parquet)")
    parser.add_argument("--spike-in-stats", type=Path, required=True, help="Input spike_in_stats.parquet")
    parser.add_argument("--output-figure", type=Path, required=True, help="Output correlation figure PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the plotting, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotCorrelationConfig(
            spike_in_stats=args.spike_in_stats,
            output_figure=args.output_figure,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
