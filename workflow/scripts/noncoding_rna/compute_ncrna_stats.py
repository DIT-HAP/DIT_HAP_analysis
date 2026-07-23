#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Non-coding RNA Stats
=====================

Stage 2a of the non-coding RNA split: read the prepared nuclear_trnas parquet
intermediate and write it straight to the stats TSV. The "stats" table IS the
per-nuclear-tRNA table itself (Systematic ID, GtRNAdb_Name, Amino_Acid,
Anticodon, tRNA_copy_number, DR/DL, mRNA abundance means) — byte-faithful to
the original single-script port's `nuclear_trnas.to_csv(config.output_stats, ...)`.
Depends only on prepare_ncrna_table's output, so it re-runs independently of
the figures rule.

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
class StatsConfig:
    """Parquet input + TSV output for the ncRNA stats table."""
    nuclear_trnas: Path
    output_stats: Path

    def validate(self) -> None:
        """Raise ValueError if the required input is missing, then ensure the output dir exists."""
        if not self.nuclear_trnas.exists():
            raise ValueError(f"Required input not found: {self.nuclear_trnas}")
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
    """Read parquet -> write the nuclear tRNA table straight to stats TSV."""
    config.validate()

    nuclear_trnas = read_parquet(config.nuclear_trnas)
    nuclear_trnas.to_csv(config.output_stats, sep="\t", index=False)

    fitted = int(nuclear_trnas["DR"].notna().sum()) if "DR" in nuclear_trnas.columns else 0
    logger.success(
        f"Non-coding RNA stats: {len(nuclear_trnas):,} nuclear tRNAs "
        f"({fitted:,} with DR), {nuclear_trnas['Anticodon'].nunique():,} distinct anticodons"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Write the nuclear tRNA table to a stats TSV")
    parser.add_argument("--nuclear-trnas", type=Path, required=True, help="Input nuclear_trnas.parquet")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output ncRNA stats TSV")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the conversion, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = StatsConfig(
            nuclear_trnas=args.nuclear_trnas,
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
