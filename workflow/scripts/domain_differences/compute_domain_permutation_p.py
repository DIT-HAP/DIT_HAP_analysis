#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute Domain Permutation P
============================

Stage 4a, step 4: the two permutation p-values behind the segment reliability
columns, as a table of their own.

They live in a separate rule because they are the only expensive part of the stage
(each gene re-runs n_permutations shuffles per segment) and because they are
reference numbers, not the analysis: the analysis is filtered by the effect floor
alone. Skipping this rule leaves ``domain_segments.tsv`` complete.

Two columns per segment, with deliberately different standing:

- ``perm_p_seg`` — this segment against the rest of its gene. ANTI-CONSERVATIVE:
  the segment was selected for exactly this contrast, so it is only comparable
  across segments produced by the same rule.
- ``perm_p_gene`` — the gene's max-over-cuts contrast against a max-over-cuts null,
  i.e. calibrated. Both are floored at ``1 / (n_permutations + 1)``. Neither is
  used as a filter: the analysis is the effect floor, and a p that only the
  selected cut can reach would put the selection back into the significance test.

Input
-----
- results/4a_domain_differences/{dataset}/_work/insertions.parquet
- results/4a_domain_differences/{dataset}/_work/segments.parquet — the index ranges

Output
------
- domain_permutation_p.tsv — ``Systematic ID, segment, perm_p_seg, perm_p_gene``

Usage
-----
    python workflow/scripts/domain_differences/compute_domain_permutation_p.py \\
        --insertions results/4a_domain_differences/HD_DIT_HAP/_work/insertions.parquet \\
        --segments results/4a_domain_differences/HD_DIT_HAP/_work/segments.parquet \\
        --n-permutations 1000 --min-seg 1 --by-site \\
        --output results/4a_domain_differences/HD_DIT_HAP/domain_permutation_p.tsv

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

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from domain_differences.core import ordered_genes  # noqa: E402
from domain_differences.segmentation import segment_permutation  # noqa: E402
from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
PERMUTATION_COLUMNS = ["Systematic ID", "segment", "perm_p_seg", "perm_p_gene"]
SEED = 0


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PermutationConfig:
    """Runtime parameters for the permutation table."""

    insertions: Path
    segments: Path
    n_permutations: int
    min_seg: int
    by_site: bool
    output: Path


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PermutationConfig) -> int:
    """Re-derive each reported gene's profile and permutation p's, then write them."""
    insertions = read_parquet(config.insertions)
    ranges = read_parquet(config.segments)
    profiles = ordered_genes(insertions, by_site=config.by_site)
    rng = np.random.default_rng(SEED)

    rows = []
    for sid, group in ranges.groupby("Systematic ID", sort=True):
        # the index ranges are positional into the profile this rule rebuilds with
        # the same by_site flag; a length mismatch means the two rules disagree
        drs = profiles[sid]["DR"].to_numpy(dtype=float)
        segments = [(int(lo), int(hi)) for lo, hi in zip(group["lo"], group["hi"])]
        if segments[-1][1] != len(drs):
            raise ValueError(
                f"{sid}: segment ranges end at {segments[-1][1]} but the profile has "
                f"{len(drs)} points — the segmentation and permutation rules disagree"
            )
        per_segment, gene_p = segment_permutation(
            drs, segments, n_perm=config.n_permutations, rng=rng, min_seg=config.min_seg
        )
        for k, p in enumerate(per_segment):
            rows.append({"Systematic ID": sid, "segment": k + 1,
                         "perm_p_seg": round(p, 4), "perm_p_gene": round(gene_p, 4)})

    table = pd.DataFrame(rows, columns=PERMUTATION_COLUMNS)
    table.to_csv(config.output, index=False, sep="\t")

    resolution = 1 / (config.n_permutations + 1)
    logger.info(f"{len(table):,} segments over {table['Systematic ID'].nunique():,} genes "
                f"({config.n_permutations:,} permutations each, floor {resolution:.4f})")
    logger.info(f"perm_p_seg  <= 0.05: {100 * table['perm_p_seg'].le(0.05).mean():.0f}%  (anti-conservative)")
    logger.info(f"perm_p_gene <= 0.05: {100 * table['perm_p_gene'].le(0.05).mean():.0f}%  (calibrated)")
    logger.info(f"wrote {config.output}")
    return 0


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Permutation p-values per segment (4a)")
    parser.add_argument("--insertions", type=Path, required=True, help="Input insertions.parquet")
    parser.add_argument("--segments", type=Path, required=True, help="Input segments.parquet")
    parser.add_argument("--n-permutations", type=int, required=True, help="Shuffles per p-value")
    parser.add_argument("--min-seg", type=int, required=True, help="Minimum points per segment")
    parser.add_argument("--by-site", action=argparse.BooleanOptionalAction, default=True,
                        help="Must match the flag the segmentation rule used")
    parser.add_argument("--output", type=Path, required=True, help="Output permutation TSV")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the permutations, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PermutationConfig(
            insertions=args.insertions.resolve(), segments=args.segments.resolve(),
            n_permutations=args.n_permutations, min_seg=args.min_seg, by_site=args.by_site,
            output=args.output.resolve(),
        )
        if config.n_permutations < 1:
            raise ValueError(f"--n-permutations must be at least 1, got {config.n_permutations}")
        config.output.parent.mkdir(parents=True, exist_ok=True)
        run(config)
    except (ValueError, FileNotFoundError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
