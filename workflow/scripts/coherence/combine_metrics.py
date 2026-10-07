#!/usr/bin/env python3

"""
Coherence Metrics Combiner (per-source tables -> one cross-source table)
========================================================================

Final stage of the coherence DAG. Concatenates the per-source
coherence_metrics.parquet tables (go_macrocomplex / go_cc / go_bp / ...) into ONE
long table so complexes, cellular components and biological processes can be
ranked, filtered and browsed side by side. Each row already carries its own
`source` column (set upstream by compute_coherence.py), so the sources stay
distinguishable after the concat.

FDR is re-derived over the UNION here, replacing the per-source q_value. A
q-value is only defined relative to a hypothesis family, and this table is a
single cross-source ranked list, so its family is its own row set. Carrying
per-source q into it meant comparing q-values computed against different family
sizes — go_bp has ~8x the terms of go_macrocomplex, so the same p maps to ~8x
different q and the smaller family won every cross-source comparison, which
deduplicate_terms.py does directly: with scope=pooled its per-cluster
representative is picked by min q over a cluster that spans sources, and the
artifact was measurable (go_macrocomplex is 8% of the terms but was 20% of the
representatives). The per-source tables themselves keep their per-source q,
which is correct for them — do NOT pool those. The per-source figures
(plot_coherence / plot_group_scatter) read those tables, not this one, so they
are unaffected by the change here.
Order matters: correct here, on the full concatenation, then let
deduplicate_terms.py carry that q through unchanged. Re-running BH on the
deduped representatives would be anti-conservative, not merely redundant: the
representative is selected as its cluster's min q, and a family of per-cluster
minima is not uniform under the null.

Input
-----
- --metrics: one or more per-source coherence_metrics.parquet paths (order is
  only cosmetic; the output is re-sorted). Each must share the same column schema.

Output
------
- --output: combined/coherence_metrics.parquet — the row-wise concatenation of
  the inputs, sorted by median_pairwise_distance_z ascending (most coherent
  first), same columns as the per-source tables except q_value, which is
  re-derived by BH over the pooled rows (see the FDR note above).

Usage
-----
    python combine_metrics.py \\
        --metrics results/3a_coherence/{dataset}/go_macrocomplex/coherence_metrics.parquet \\
                  results/3a_coherence/{dataset}/go_cc/coherence_metrics.parquet \\
                  results/3a_coherence/{dataset}/go_bp/coherence_metrics.parquet \\
        --output results/3a_coherence/{dataset}/combined/coherence_metrics.parquet

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-09-20
Version:  1.1.0
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

# 3. Third-party Imports
from loguru import logger
from scipy.stats import false_discovery_control

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from io_table import read_parquet, write_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# The statistic the pooled FDR correction runs on: compute_coherence.py's primary
# coherence method, the same p the per-source q_value was derived from. Spelled
# out here rather than imported — that module is a CLI script, and this string is
# the cross-script schema contract (the dedup stage reads q_value, not this).
_PRIMARY_P_COLUMN = "median_pairwise_distance_p"


# =============================================================================
# CORE LOGIC
# =============================================================================
def combine(metrics_paths: list[Path]) -> pd.DataFrame:
    """Concat per-source metrics, re-correct q over the union, sort by primary z ascending."""
    # Empty per-source tables (a source where no group passed the size filter) are
    # tolerated and contribute no rows. If every input is empty the result is an
    # empty frame with no columns, which the Parquet writer round-trips as a 0x0
    # table — downstream consumers already handle the empty-table case.
    frames = []
    for path in metrics_paths:
        frame = read_parquet(path)
        if frame.empty:
            logger.warning(f"empty metrics table (no groups passed the filter): {path}")
            continue
        frames.append(frame)
    if not frames:
        logger.warning("all per-source metrics tables were empty; writing an empty combined table")
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    # The pooled family's own correction, replacing every per-source q. Indexing the
    # primary method's p column directly (no has-column guard) is deliberate: a
    # missing column is a schema break upstream, and raising here beats silently
    # leaving the per-source q in place under the combined table's name.
    combined["q_value"] = false_discovery_control(
        combined[_PRIMARY_P_COLUMN].to_numpy(), method="bh"
    )
    if "median_pairwise_distance_z" in combined.columns:
        combined = combined.sort_values("median_pairwise_distance_z").reset_index(drop=True)
    return combined


@logger.catch(reraise=True)
def run(metrics_paths: list[Path], output: Path) -> None:
    """Combine the per-source metrics tables and write the unified table."""
    output.parent.mkdir(parents=True, exist_ok=True)
    combined = combine(metrics_paths)
    write_parquet(combined, output)
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
    parser.add_argument("--metrics", type=Path, nargs="+", required=True, help="Per-source coherence_metrics.parquet paths")
    parser.add_argument("--output", type=Path, required=True, help="Output combined/coherence_metrics.parquet")
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
