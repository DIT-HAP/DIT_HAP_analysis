#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
UTR Insertion Classification
============================

Stage 2 of the UTR split (see workflow/rules/utr.smk's header comment): reads
the three parquet intermediates written by prepare_utr_data.py, classifies
intergenic insertions within UTR_DISTANCE_THRESHOLD (400 bp) of a flanking
gene boundary as 5'UTR or 3'UTR insertions (strand-aware), and writes the
per-insertion UTR stats TSV with um_ratio (insertion DR / gene DR) and A_ratio
(insertion A / gene A). Ported from the deterministic section of
DIT_HAP_pipeline/workflow/notebooks/upstream_and_downstream_analysis.ipynb.

The downstream human-review / plotting notebook lives at
notebooks/domain_analysis/review_utr_insertions.ipynb (Task 9) — this script
only produces the deterministic per-insertion table.

Input
-----
- fitting_results.parquet: insertion-level A/DR, indexed by [Chr, Coordinate,
  Strand, Target] (from prepare_utr_data.py).
- annotations.parquet: insertion-level annotations (same MultiIndex). For
  intergenic rows, Name / Systematic ID / Strand_Interval are pipe-separated
  "left|right" pairs describing the two genes flanking the intergenic
  interval, and Distance_to_region_start / _end are the distances to the LEFT
  gene's 3' end and the RIGHT gene's 5' end respectively.
- gene_result.parquet: gene-level fitting stats (Systematic ID, Name,
  DeletionLibrary_essentiality, A, DR, ...). Joined to insertions on the short
  gene Name (the annotations' Name field carries short names).

Output
------
- utr_insertion_stats.tsv: one row per UTR insertion (an intergenic insertion
  assigned to a flanking gene), with Parental_gene, UTR_type (5UTR/3UTR),
  Insertion_direction, Distance_to_gene_boundary, insertion/gene A + DR, and
  the derived um_ratio + A_ratio.

Usage
-----
    python classify_utr_insertions.py \\
        --fitting-results results/utr/{dataset}/_work/fitting_results.parquet \\
        --annotations results/utr/{dataset}/_work/annotations.parquet \\
        --gene-result results/utr/{dataset}/_work/gene_result.parquet \\
        --distance-threshold 400 \\
        --output-stats results/utr/{dataset}/utr_insertion_stats.tsv

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-22
Version:  2.0.0
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
from workflow.src.utr.core import UTR_DISTANCE_THRESHOLD, classify_utr_insertions  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class UTRConfig:
    """Inputs, outputs, and threshold for the UTR insertion classification."""
    fitting_results: Path
    annotations: Path
    gene_result: Path
    output_stats: Path
    distance_threshold: int = UTR_DISTANCE_THRESHOLD

    def validate(self) -> None:
        """Raise ValueError if any input is missing or the threshold is invalid, then ensure the output dir exists."""
        for path in [self.fitting_results, self.annotations, self.gene_result]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.distance_threshold <= 0:
            raise ValueError(f"distance_threshold must be positive, got {self.distance_threshold}")
        self.output_stats.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: UTRConfig) -> None:
    """Load parquets -> classify UTR insertions -> save per-insertion TSV."""
    config.validate()

    fitting_results = read_parquet(config.fitting_results)
    annotations = read_parquet(config.annotations)
    gene_result = read_parquet(config.gene_result)

    stats = classify_utr_insertions(fitting_results, annotations, gene_result, config.distance_threshold)
    stats.to_csv(config.output_stats, sep="\t", index=False)

    n5 = int((stats["UTR_type"] == "5UTR").sum())
    n3 = int((stats["UTR_type"] == "3UTR").sum())
    logger.success(
        f"UTR insertions classified: {len(stats):,} total "
        f"({n5:,} 5'UTR, {n3:,} 3'UTR) -> {config.output_stats}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Classify intergenic insertions near gene boundaries as 5'/3' UTR insertions")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Insertion-level fitting_results.parquet")
    parser.add_argument("--annotations", type=Path, required=True, help="Insertion-level annotations.parquet")
    parser.add_argument("--gene-result", type=Path, required=True, help="Gene-level gene_result.parquet")
    parser.add_argument("--distance-threshold", type=int, default=UTR_DISTANCE_THRESHOLD, help="UTR distance threshold in bp (default: 400)")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output per-insertion UTR stats TSV")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = UTRConfig(
            fitting_results=args.fitting_results,
            annotations=args.annotations,
            gene_result=args.gene_result,
            output_stats=args.output_stats,
            distance_threshold=args.distance_threshold,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
