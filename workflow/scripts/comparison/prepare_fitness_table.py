#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare Fitness Table
======================

Stage 1 of the comparison split: merge the PomBase-derived protein-features
table (the spine, carrying the other large-scale studies' fitness/depletion
columns) with the DIT-HAP and gRNA fitness metrics from the gene annotation
reference parquet (1c_annotate), producing the fitness_table.parquet
intermediate consumed by the stats / figures rules. Integration-density
columns are clipped at --clip-upper on the way out.

Input
-----
- pombe_coding_gene_protein_features.tsv (key gene_systematic_id; Barseq,
  integration density, ipkm/uipkm, colony size, growth-rate columns).
- gene_annotation_reference.protein.parquet (index gene_systematic_id;
  HD_DIT_HAP_DR and gRNA_DR metric columns).

Output
------
- fitness_table.parquet: features spine left-joined with the two DR metrics,
  density columns clipped.

Usage
-----
    python workflow/scripts/comparison/prepare_fitness_table.py \\
        --protein-features results/1b_features/2026-06-01/pombe_coding_gene_protein_features.tsv \\
        --annotation-reference results/1c_annotation/2026-06-01/2026-08-11/gene_annotation_reference.protein.parquet \\
        --output-fitness-table results/6a_comparison/HD_DIT_HAP/_work/fitness_table.parquet

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-10-08
Version:  3.0.0
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

# 4. Local Imports (relative path resolution)
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from io_table import read_parquet, write_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402
from comparison.core import (  # noqa: E402
    CLIP_UPPER,
    build_fitness_table,
)


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True, slots=True)
class PrepareConfig:
    """Feature TSV + annotation reference parquet inputs and the parquet output."""
    protein_features: Path
    annotation_reference: Path
    clip_upper: float
    output_fitness_table: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dir exists."""
        for path in [self.protein_features, self.annotation_reference]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_fitness_table.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PrepareConfig) -> None:
    """Load -> merge -> clip -> write the fitness table parquet intermediate."""
    config.validate()

    protein_features = pd.read_csv(config.protein_features, sep="\t")
    annotation_reference = read_parquet(config.annotation_reference)

    fitness_table = build_fitness_table(
        annotation_reference, protein_features, clip_upper=config.clip_upper
    )
    write_parquet(fitness_table, config.output_fitness_table)

    logger.success(f"Prepared fitness table: {len(fitness_table):,} genes")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare the merged fitness table parquet intermediate")
    parser.add_argument("--protein-features", type=Path, required=True, help="pombe_coding_gene_protein_features.tsv")
    parser.add_argument("--annotation-reference", type=Path, required=True, help="gene_annotation_reference.protein.parquet")
    parser.add_argument("--clip-upper", type=float, default=CLIP_UPPER, help="Upper cap for integration-density columns")
    parser.add_argument("--output-fitness-table", type=Path, required=True, help="Output fitness_table.parquet")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the preparation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PrepareConfig(
            protein_features=args.protein_features,
            annotation_reference=args.annotation_reference,
            clip_upper=args.clip_upper,
            output_fitness_table=args.output_fitness_table,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
