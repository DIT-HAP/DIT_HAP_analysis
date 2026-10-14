#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare Domain Insertions
=========================

Stage 4a, step 1: turn the upstream release's three tables into one positioned
insertion profile. Every in-gene insertion of every gene gets a row carrying its
position along the CDS in both coordinate systems the stage uses — ``frac`` for
the prototype's legacy columns and ``Residue_affected`` for everything the
gallery draws — plus the insertion site key (``Coordinate``) the by-site collapse
groups on and the two directions' DR values.

The expensive part of the stage is the annotation dedup + reindex; doing it once
here is what lets the segmentation rule stay a pure function of this parquet.

Input
-----
- release ``insertion_level/fitting_results.tsv`` — the insertion universe.
- release ``insertion_level/annotations.tsv.gz`` — in-gene call, CDS distances.
- release ``gene_level/fitting_results.tsv`` — gene names.

Output
------
- ``_work/insertions.parquet`` — one row per in-gene insertion, sorted by
  (Systematic ID, frac). Columns: ``Systematic ID, name, Chr, Coordinate, Strand,
  Insertion_direction, Residue_affected, Residue_frame, frac, DR, DL``.

Usage
-----
    python workflow/scripts/domain_differences/prepare_domain_insertions.py \\
        --fitting-results <release>/insertion_level/fitting_results.tsv \\
        --annotations <release>/insertion_level/annotations.tsv.gz \\
        --gene-level <release>/gene_level/fitting_results.tsv \\
        --output results/4a_domain_differences/HD_DIT_HAP/_work/insertions.parquet

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-10-14
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
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from domain_differences.core import read_insertions  # noqa: E402
from io_table import write_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PrepareConfig:
    """Runtime parameters for the insertion-profile preparation."""

    fitting_results: Path
    annotations: Path
    gene_level: Path
    output: Path


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PrepareConfig) -> int:
    """Build the positioned insertion profile and write it as parquet."""
    insertions = read_insertions(config.fitting_results, config.annotations, config.gene_level)
    write_parquet(insertions, config.output)

    genes = insertions["Systematic ID"].nunique()
    sites = insertions.groupby(["Systematic ID", "Coordinate"], observed=True).ngroups
    logger.info(
        f"{len(insertions):,} in-gene insertions over {genes:,} genes "
        f"({sites:,} insertion sites, {len(insertions) / sites:.2f} insertions per site)"
    )
    logger.info(f"wrote {config.output}")
    return 0


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare the positioned insertion profile for 4a")
    parser.add_argument("--fitting-results", type=Path, required=True,
                        help="Release insertion_level/fitting_results.tsv")
    parser.add_argument("--annotations", type=Path, required=True,
                        help="Release insertion_level/annotations.tsv.gz")
    parser.add_argument("--gene-level", type=Path, required=True,
                        help="Release gene_level/fitting_results.tsv")
    parser.add_argument("--output", type=Path, required=True, help="Output parquet path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the preparation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PrepareConfig(
            fitting_results=args.fitting_results.resolve(),
            annotations=args.annotations.resolve(),
            gene_level=args.gene_level.resolve(),
            output=args.output.resolve(),
        )
        config.output.parent.mkdir(parents=True, exist_ok=True)
        run(config)
    except (ValueError, FileNotFoundError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
