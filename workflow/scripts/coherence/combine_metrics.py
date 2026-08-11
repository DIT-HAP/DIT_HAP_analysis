#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coherence Metrics Combiner (per-source tables -> one cross-source table)
========================================================================

Final stage of the coherence DAG. Concatenates the per-source
coherence_metrics.tsv tables (go_macrocomplex / go_cc / go_bp / ...) into ONE
long table so complexes, cellular components and biological processes can be
ranked, filtered and browsed side by side. Each row already carries its own
`source` column (set upstream by compute_coherence.py), so the sources stay
distinguishable after the concat.

No re-correction happens here: p_fdr is computed per source in
compute_coherence.py (each source is its own hypothesis family), and pooling the
already-corrected rows would conflate families. This script only stacks and
sorts; the per-source p_fdr is carried through unchanged.

Input
-----
- --metrics: one or more per-source coherence_metrics.tsv paths (order is only
  cosmetic; the output is re-sorted). Each must share the same column schema.

Output
------
- --output: coherence_metrics_combined.tsv — the row-wise concatenation of the
  inputs, sorted by z_score ascending (most coherent first), same columns as the
  per-source tables.

Usage
-----
    python combine_metrics.py \\
        --metrics results/coherence/{dataset}/go_macrocomplex/coherence_metrics.tsv \\
                  results/coherence/{dataset}/go_cc/coherence_metrics.tsv \\
                  results/coherence/{dataset}/go_bp/coherence_metrics.tsv \\
        --output results/coherence/{dataset}/coherence_metrics_combined.tsv

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-23
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd
from pandas.errors import EmptyDataError

# 3. Third-party Imports
from loguru import logger


# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
def combine(metrics_paths: list[Path]) -> pd.DataFrame:
    """Row-concatenate per-source metrics tables, sorted by z_score ascending.

    Empty per-source tables (a source where no group passed the size filter)
    are tolerated and contribute no rows. If every input is empty the result is
    an empty frame with no columns, which to_csv still writes as a valid header-
    less TSV — downstream consumers already handle the empty-table case.
    """
    frames = []
    for path in metrics_paths:
        try:
            frames.append(pd.read_csv(path, sep="\t"))
        except EmptyDataError:
            # A source where no group passed the size filter writes a truly
            # empty TSV (compute_coherence.py's to_csv on a 0-column frame emits
            # no header line). Skip it — it contributes no rows.
            logger.warning(f"empty metrics table (no groups passed the filter): {path}")
    non_empty = [f for f in frames if not f.empty]
    if not non_empty:
        logger.warning("all per-source metrics tables were empty; writing an empty combined table")
        return pd.DataFrame()
    combined = pd.concat(non_empty, ignore_index=True)
    if "z_score" in combined.columns:
        combined = combined.sort_values("z_score").reset_index(drop=True)
    return combined


@logger.catch(reraise=True)
def run(metrics_paths: list[Path], output: Path) -> None:
    """Combine the per-source metrics tables and write the unified table."""
    output.parent.mkdir(parents=True, exist_ok=True)
    combined = combine(metrics_paths)
    combined.to_csv(output, sep="\t", index=False)
    by_source = (
        combined["source"].value_counts().to_dict() if "source" in combined.columns else {}
    )
    logger.success(f"combined {len(metrics_paths)} sources -> {len(combined):,} groups "
                   f"({by_source}) -> {output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Combine per-source coherence metrics into one cross-source table")
    parser.add_argument("--metrics", type=Path, nargs="+", required=True, help="Per-source coherence_metrics.tsv paths")
    parser.add_argument("--output", type=Path, required=True, help="Output coherence_metrics_combined.tsv")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: parse args, combine the per-source tables, write the unified table."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        run(args.metrics, args.output)
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
