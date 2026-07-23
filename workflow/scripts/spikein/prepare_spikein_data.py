#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare Spike-In Data
======================

Stage 1 of the spike-in split: load the filtered raw-reads insertion table,
extract the known spike-in sites, and assign each sample its dilution ratio
by read-count rank, writing the long-form stats table as a parquet
intermediate consumed by the compute-stats / plot-correlation rules.

Input
-----
- Filtered raw-reads insertion table (Spikein's pre-release results/13_filtered/
  raw_reads.filtered.tsv — release/ never packages this file, see spikein.smk),
  indexed by [Chr, Coordinate, Strand] (a 4th Target level, if present, is
  dropped), columned by [Sample, Timepoint] (one Timepoint per dilution point).

Output
------
- spike_in_stats.parquet: long-form per-site/per-sample table (Reads, Ratio,
  Relative_Read_Ratio, Relative_Dilution_Ratio).

Usage
-----
    python prepare_spikein_data.py \\
        --raw-reads .../Spikein/results/13_filtered/raw_reads.filtered.tsv \\
        --output-spike-in-stats results/spikein/_work/spike_in_stats.parquet \\
        --spike-in-sites-json '{"DY215": {"chr": "I", "coord": 3749394, "strand": "-"}, ...}'

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-22
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io import write_parquet  # noqa: E402
from workflow.src.spikein.core import DEFAULT_SPIKE_IN_SITES, build_spike_in_stats  # noqa: E402


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PrepareConfig:
    """Input raw-reads table, spike-in site coordinates, and the parquet output."""
    raw_reads: Path
    output_spike_in_stats: Path
    spike_in_sites: dict[str, dict]

    def validate(self) -> None:
        """Raise ValueError if the input is missing, then ensure the output dir exists."""
        if not self.raw_reads.exists():
            raise ValueError(f"Required input not found: {self.raw_reads}")
        self.output_spike_in_stats.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PrepareConfig) -> None:
    """Load raw reads -> extract spike-in sites -> assign ratios -> write parquet."""
    config.validate()

    raw_reads = pd.read_csv(config.raw_reads, sep="\t", header=[0, 1], index_col=[0, 1, 2, 3])
    spike_in_stats = build_spike_in_stats(raw_reads, config.spike_in_sites)

    write_parquet(spike_in_stats, config.output_spike_in_stats)

    logger.success(f"Prepared spike-in stats: {len(spike_in_stats):,} site x sample rows")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare spike-in stats parquet intermediate")
    parser.add_argument("--raw-reads", type=Path, required=True, help="Filtered raw-reads insertion table (tsv)")
    parser.add_argument("--output-spike-in-stats", type=Path, required=True, help="Output spike_in_stats.parquet")
    parser.add_argument(
        "--spike-in-sites-json", type=str, default=None,
        help="JSON dict of {strain: {chr, coord, strand}} (default: the 5 hardcoded DY sites)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the preparation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        spike_in_sites = (
            json.loads(args.spike_in_sites_json) if args.spike_in_sites_json else DEFAULT_SPIKE_IN_SITES
        )
        config = PrepareConfig(
            raw_reads=args.raw_reads,
            output_spike_in_stats=args.output_spike_in_stats,
            spike_in_sites=spike_in_sites,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
