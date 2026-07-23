#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare Non-coding RNA Tables
==============================

Stage 1 of the non-coding RNA split: load the non-coding-gene fitting results,
ncRNA genome-region bed, GtRNAdb tRNA annotations, and Marguerat 2012 mRNA
abundance, then merge them into two parquet intermediates consumed by the
stats / figures rules:

- combined.parquet: the full annotated ncRNA table (bed + GtRNAdb + fitting +
  abundance), needed by plot_feature_type_donut.
- nuclear_trnas.parquet: nuclear tRNAs only, annotated with amino
  acid/anticodon/copy number, needed by both the stats TSV and
  plot_trna_summary.

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
from loguru import logger

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io import write_parquet  # noqa: E402
from workflow.src.noncoding_rna.core import (  # noqa: E402
    build_ncrna_table,
    load_gtrnadb,
    load_marguerat_abundance,
    load_ncrna_fitting,
    select_nuclear_tRNAs,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PrepareConfig:
    """Inputs and parquet outputs for the ncRNA table preparation."""
    ncrna_fitting: Path
    ncrna_bed: Path
    gtrnadb_bed: Path
    marguerat_excel: Path
    output_combined: Path
    output_nuclear_trnas: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.ncrna_fitting, self.ncrna_bed, self.gtrnadb_bed, self.marguerat_excel]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_combined, self.output_nuclear_trnas]:
            out.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PrepareConfig) -> None:
    """Load -> merge -> select nuclear tRNAs -> write the two parquet intermediates."""
    config.validate()

    ncrna_fitting = load_ncrna_fitting(config.ncrna_fitting)
    ncrna_bed = pd.read_csv(config.ncrna_bed, sep="\t")
    gtrnadb = load_gtrnadb(config.gtrnadb_bed)
    marguerat_means = load_marguerat_abundance(config.marguerat_excel)

    combined = build_ncrna_table(ncrna_fitting, ncrna_bed, gtrnadb, marguerat_means)
    nuclear_trnas = select_nuclear_tRNAs(combined)

    write_parquet(combined, config.output_combined)
    write_parquet(nuclear_trnas, config.output_nuclear_trnas)

    logger.success(
        f"Prepared ncRNA tables: {len(combined):,} annotated ncRNA genes, "
        f"{len(nuclear_trnas):,} nuclear tRNAs"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare non-coding RNA parquet intermediates")
    parser.add_argument("--ncrna-fitting", type=Path, required=True, help="Non-coding-gene fitting_results.tsv")
    parser.add_argument("--ncrna-bed", type=Path, required=True, help="ncRNA genome-region bed")
    parser.add_argument("--gtrnadb-bed", type=Path, required=True, help="GtRNAdb tRNA bed (schiPomb_972H-tRNAs.bed)")
    parser.add_argument("--marguerat-excel", type=Path, required=True, help="Marguerat 2012 abundance xlsx")
    parser.add_argument("--output-combined", type=Path, required=True, help="Output combined.parquet")
    parser.add_argument("--output-nuclear-trnas", type=Path, required=True, help="Output nuclear_trnas.parquet")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the preparation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PrepareConfig(
            ncrna_fitting=args.ncrna_fitting,
            ncrna_bed=args.ncrna_bed,
            gtrnadb_bed=args.gtrnadb_bed,
            marguerat_excel=args.marguerat_excel,
            output_combined=args.output_combined,
            output_nuclear_trnas=args.output_nuclear_trnas,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
