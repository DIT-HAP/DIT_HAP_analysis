#!/usr/bin/env python3

"""
Compute LFC Coverage
====================

Per-dataset half of the LFC-coverage comparison: flag every protein-coding gene
as LFC-covered (an insertion inside it survived the release's hard filter and
replicate merge, i.e. it has a row in insertion_level/LFC.tsv and passes
IN_GENE_FILTER) and DR-covered (DR not NaN), then summarize both by FYPO
viability class. The per-gene flags are kept as a parquet so the statistics can
be audited without re-reading the release.

Input
-----
- {release_dir}/insertion_level/LFC.tsv
- {release_dir}/insertion_level/annotations.tsv.gz
- results/2a_coverage/{dataset}/_work/gene_result.parquet

Output
------
- results/2a_coverage/{dataset}/lfc_coverage_stats.tsv
- results/2a_coverage/{dataset}/_work/lfc_gene_coverage.parquet

Usage
-----
    python workflow/scripts/coverage/compute_lfc_coverage.py \\
        --lfc .../insertion_level/LFC.tsv \\
        --annotations .../insertion_level/annotations.tsv.gz \\
        --gene-result results/2a_coverage/HD_DIT_HAP/_work/gene_result.parquet \\
        --dataset HD_DIT_HAP \\
        --output-stats results/2a_coverage/HD_DIT_HAP/lfc_coverage_stats.tsv \\
        --output-gene-table results/2a_coverage/HD_DIT_HAP/_work/lfc_gene_coverage.parquet

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-07
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

# 2. Third-party Imports
from loguru import logger
import pandas as pd

# 3. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from coverage.core import build_lfc_gene_table, summarize_lfc_coverage
from io_table import read_parquet, write_parquet
from logging_setup import setup_logger

# =============================================================================
# CONSTANTS
# =============================================================================
# gene_result's boundary contract: the gene universe, the class column, and the
# DR column the reference statistic is read from.
REQUIRED_GENE_RESULT_COLUMNS = ("Systematic ID", "FYPOviability", "DR")


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class ComputeLfcCoverageConfig:
    """Release LFC tables, the prepared gene universe, and the two outputs."""

    lfc: Path
    annotations: Path
    gene_result: Path
    dataset: str
    output_stats: Path
    output_gene_table: Path

    def validate(self) -> None:
        """Raise ValueError if an input is missing, then ensure output dirs exist."""
        for path in [self.lfc, self.annotations, self.gene_result]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_stats.parent.mkdir(parents=True, exist_ok=True)
        self.output_gene_table.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
def run(config: ComputeLfcCoverageConfig) -> None:
    """Flag per-gene LFC/DR coverage, summarize by FYPO class, write both tables."""
    config.validate()
    gene_result = read_parquet(config.gene_result)
    # Narrow read_parquet's honest DataFrame | Series return for the type checker;
    # this stage's input is always a table, so a Series is a wrong-file error.
    if not isinstance(gene_result, pd.DataFrame):
        raise TypeError(f"{config.gene_result} must be a table, not a Series")

    missing_columns = set(REQUIRED_GENE_RESULT_COLUMNS) - set(gene_result.columns)
    if missing_columns:
        raise ValueError(
            f"{config.gene_result} is missing required columns: {sorted(missing_columns)}"
        )

    gene_table = build_lfc_gene_table(config.lfc, config.annotations, gene_result)
    stats = summarize_lfc_coverage(gene_table)
    stats.insert(0, "dataset", config.dataset)

    stats.to_csv(config.output_stats, sep="\t", index=False)
    write_parquet(gene_table, config.output_gene_table)

    total = stats[stats["fypo_class"] == "total"].iloc[0]
    n_without_fypo = int(gene_table["FYPOviability"].isna().sum())
    logger.success(
        f"{config.dataset}: LFC-covered {total['n_lfc']:,}/{total['n_total']:,} "
        f"({total['pct_lfc']:.1%}); DR-covered {total['n_dr']:,} "
        f"({total['pct_dr']:.1%}); {n_without_fypo:,} genes without FYPO annotation"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Compute per-dataset LFC gene coverage")
    parser.add_argument("--lfc", type=Path, required=True, help="Insertion-level LFC.tsv")
    parser.add_argument(
        "--annotations", type=Path, required=True, help="Insertion-level annotations.tsv.gz"
    )
    parser.add_argument(
        "--gene-result", type=Path, required=True, help="Prepared gene_result.parquet"
    )
    parser.add_argument("--dataset", required=True, help="Dataset name stamped on every stats row")
    parser.add_argument(
        "--output-stats", type=Path, required=True, help="Output lfc_coverage_stats.tsv"
    )
    parser.add_argument(
        "--output-gene-table", type=Path, required=True, help="Output per-gene coverage parquet"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging"
    )
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the computation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = ComputeLfcCoverageConfig(
            lfc=args.lfc,
            annotations=args.annotations,
            gene_result=args.gene_result,
            dataset=args.dataset,
            output_stats=args.output_stats,
            output_gene_table=args.output_gene_table,
        )
        run(config)
    except Exception as e:
        logger.exception(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
