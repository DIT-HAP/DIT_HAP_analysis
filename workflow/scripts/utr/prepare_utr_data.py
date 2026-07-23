#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare UTR Data
================

Stage 1 of the UTR split (see workflow/rules/utr.smk's header comment): load
the insertion-level fitting_results + annotations and the gene-level
fitting_results, normalizing legacy um/lam -> DR/DL columns, then write three
parquet intermediates consumed by the classification rule:

- fitting_results.parquet: insertion-level per-insertion A/DR, indexed by
  [Chr, Coordinate, Strand, Target].
- annotations.parquet: insertion-level annotations (same MultiIndex),
  including the pipe-separated two-gene intergenic-interval columns.
- gene_result.parquet: gene-level fitting stats (Systematic ID, Name,
  DeletionLibrary_essentiality, A, DR, ...).

Author:   Yusheng Yang (guidance) + Claude (implementation)
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
from workflow.src.io import write_parquet  # noqa: E402
from workflow.src.utr.core import load_gene_level, load_insertion_level  # noqa: E402


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PrepareConfig:
    """Inputs and parquet outputs for the UTR data preparation."""
    fitting_results: Path
    annotations: Path
    gene_level: Path
    output_fitting_results: Path
    output_annotations: Path
    output_gene_result: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.fitting_results, self.annotations, self.gene_level]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_fitting_results, self.output_annotations, self.output_gene_result]:
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
    """Load -> normalize -> write the three parquet intermediates."""
    config.validate()

    fitting_results, annotations = load_insertion_level(config.fitting_results, config.annotations)
    gene_result = load_gene_level(config.gene_level)

    write_parquet(fitting_results, config.output_fitting_results)
    write_parquet(annotations, config.output_annotations)
    write_parquet(gene_result, config.output_gene_result)

    logger.success(
        f"Prepared UTR intermediates: {len(fitting_results):,} insertion fitting rows, "
        f"{len(annotations):,} annotation rows, {len(gene_result):,} genes"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare UTR parquet intermediates")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Insertion-level fitting_results.tsv")
    parser.add_argument("--annotations", type=Path, required=True, help="Insertion-level annotations.tsv(.gz)")
    parser.add_argument("--gene-level", type=Path, required=True, help="Gene-level fitting_results.tsv")
    parser.add_argument("--output-fitting-results", type=Path, required=True, help="Output fitting_results.parquet")
    parser.add_argument("--output-annotations", type=Path, required=True, help="Output annotations.parquet")
    parser.add_argument("--output-gene-result", type=Path, required=True, help="Output gene_result.parquet")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the preparation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PrepareConfig(
            fitting_results=args.fitting_results,
            annotations=args.annotations,
            gene_level=args.gene_level,
            output_fitting_results=args.output_fitting_results,
            output_annotations=args.output_annotations,
            output_gene_result=args.output_gene_result,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
