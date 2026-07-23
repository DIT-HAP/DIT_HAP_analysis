#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute Spike-In Stats
========================

Stage 2a of the spike-in split: read the prepared spike_in_stats parquet
intermediate and emit it as the long-form stats TSV. The parquet written by
prepare_spikein_data.py already IS the final stats table (site/sample-level
Reads, Ratio, Relative_Read_Ratio, Relative_Dilution_Ratio); this stage's job
is the format conversion to the release-facing TSV. Depends only on
prepare_spikein_data's output, so it re-runs independently of the
correlation-plot rule.

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


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class ComputeStatsConfig:
    """Parquet input + TSV output for the spike-in stats conversion."""
    spike_in_stats: Path
    output_stats: Path

    def validate(self) -> None:
        """Raise ValueError if the input is missing, then ensure the output dir exists."""
        if not self.spike_in_stats.exists():
            raise ValueError(f"Required input not found: {self.spike_in_stats}")
        self.output_stats.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: ComputeStatsConfig) -> None:
    """Read parquet -> write the long-form stats TSV."""
    config.validate()

    spike_in_stats = read_parquet(config.spike_in_stats)
    spike_in_stats.to_csv(config.output_stats, sep="\t")

    logger.success(f"Spike-in stats: {len(spike_in_stats):,} site x sample rows")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Spike-in stats TSV (from the prepared parquet)")
    parser.add_argument("--spike-in-stats", type=Path, required=True, help="Input spike_in_stats.parquet")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output spike_in_stats.tsv")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the conversion, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = ComputeStatsConfig(
            spike_in_stats=args.spike_in_stats,
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
