#!/usr/bin/env python3

"""
Global Coherence Null Cache (one permutation null per group size)
=================================================================

The permutation test in compute_coherence.py draws `n_scored_members` genes from
the DR-filtered background and recomputes the median/mean pairwise distance of
each draw. Those draws depend only on the group SIZE, the background cloud and
the seed — every group of the same size sees the same null. Recomputing it per
group repeated the same work thousands of times; at 100,000 permutations that is
the difference between hours and minutes.

This script computes the null once for every group size that any registered
source will score, across all sources, and writes one `.npz` cache. The
per-source compute rules then read the cache instead of drawing. The cache is
exact: it uses the same `default_rng(random_state)` + `choice` sequence as
`coherence.metrics._permutation_nulls`.

Input
-----
- fitting_results.tsv: the upstream per-gene DR/DL table (defines the background).
- group_annotation_long.tsv: one prepared long table per source; only the group
  membership and `n_annotated_members` are read, to reproduce the same filter
  compute_coherence.py applies.

Output
------
- --output: an .npz with `meta` (JSON) + `median_<n>` / `mean_<n>` arrays of
  length n_permutations for every distinct n_scored_members.

Usage
-----
    python build_null_cache.py \\
        --fitting-results .../fitting_results.tsv \\
        --annotation results/3a_coherence/{dataset}/go_bp/group_annotation_long.tsv \\
                     results/3a_coherence/{dataset}/go_cc/group_annotation_long.tsv \\
        --min-size 3 --max-size 300 --max-term-genes 350 --dr-threshold -0.3 \\
        --n-permutations 100000 --random-state 42 --jobs 8 \\
        --output results/3a_coherence/{dataset}/_null_cache/n100000_seed42.npz

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-10-20
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
import time

from loguru import logger

# 2. Third-party Imports
import pandas as pd

# 3. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from coherence.groups import build_groups
from coherence.io import load_fitting_results, load_long_table
from coherence.null_cache import (
    METHOD_KEYS,
    background_fingerprint,
    build_null_cache,
    save_null_cache,
)
from logging_setup import setup_logger


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class CacheConfig:
    """Inputs and parameters for the global null cache."""

    fitting_results: Path
    annotations: list[Path]
    output: Path
    min_size: int = 3
    max_size: int = 300
    max_term_genes: int = 500
    dr_threshold: float = -0.3
    n_permutations: int = 1000
    random_state: int = 42
    jobs: int = 1

    def validate(self) -> None:
        """Raise ValueError on missing inputs or invalid parameters."""
        required = [self.fitting_results, *self.annotations]
        for path in required:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if not self.annotations:
            raise ValueError("at least one --annotation table is required")
        if self.min_size < 2:
            raise ValueError(f"min_size must be >= 2: {self.min_size}")
        if self.max_size < self.min_size:
            raise ValueError(f"max_size ({self.max_size}) < min_size ({self.min_size})")
        if self.n_permutations < 1:
            raise ValueError(f"n_permutations must be >= 1: {self.n_permutations}")
        if self.jobs < 1:
            raise ValueError(f"jobs must be >= 1: {self.jobs}")


# =============================================================================
# CORE LOGIC
# =============================================================================
def collect_group_sizes(config: CacheConfig, background: pd.DataFrame) -> list[int]:
    """Every n_scored_members any source will score, under the compute filters."""
    sizes: set[int] = set()
    for path in config.annotations:
        long_table = load_long_table(path)
        groups = build_groups(
            background,
            long_table,
            config.min_size,
            config.max_size,
            config.max_term_genes,
        )
        sizes.update(len(grp) for grp in groups.values())
    return sorted(sizes)


@logger.catch(reraise=True)
def run(config: CacheConfig) -> None:
    """Load -> enumerate group sizes -> build the global cache -> write .npz."""
    config.validate()
    fitting = load_fitting_results(config.fitting_results)
    background = fitting[fitting["DR"] < config.dr_threshold].reset_index(drop=True)
    background_points = background[["norm_DR", "norm_DL"]].to_numpy(dtype=float)
    logger.info(
        f"background: {len(background):,} genes with DR < {config.dr_threshold}; "
        f"{len(config.annotations)} annotation table(s)"
    )

    n_values = collect_group_sizes(config, background)
    if not n_values:
        raise ValueError("no group sizes survived the filters; nothing to cache")
    logger.info(
        f"group sizes to cache: {len(n_values):,} distinct n "
        f"({n_values[0]}..{n_values[-1]}), {config.n_permutations:,} permutations each, "
        f"{config.jobs} job(s)"
    )

    started = time.perf_counter()
    nulls = build_null_cache(
        background_points, n_values, config.n_permutations, config.random_state, config.jobs
    )
    elapsed = time.perf_counter() - started

    meta = {
        "n_permutations": config.n_permutations,
        "random_state": config.random_state,
        "background_n": int(background_points.shape[0]),
        "background_fingerprint": background_fingerprint(background_points),
        "dr_threshold": config.dr_threshold,
        "min_size": config.min_size,
        "max_size": config.max_size,
        "max_term_genes": config.max_term_genes,
        "methods": list(METHOD_KEYS),
        "n_values": n_values,
        "version": 1,
    }
    save_null_cache(config.output, nulls, meta)
    logger.success(
        f"wrote {config.output} ({config.output.stat().st_size / 1024 / 1024:.1f} MB): "
        f"{len(n_values)} sizes x {config.n_permutations:,} draws in {elapsed / 60:.1f} min"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Build the global coherence null cache")
    parser.add_argument(
        "--fitting-results",
        type=Path,
        required=True,
        help="Upstream fitting_results.tsv (defines the background)",
    )
    parser.add_argument(
        "--annotation",
        type=Path,
        nargs="+",
        required=True,
        help="Prepared group_annotation_long.tsv, one per source",
    )
    parser.add_argument("--min-size", type=int, default=3)
    parser.add_argument("--max-size", type=int, default=300)
    parser.add_argument("--max-term-genes", type=int, default=500)
    parser.add_argument("--dr-threshold", type=float, default=-0.3)
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--jobs", type=int, default=1, help="Processes used to fill the cache (one per group size)"
    )
    parser.add_argument("--output", type=Path, required=True, help="Output .npz null cache")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def main() -> int:
    """Build config, run, report."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = CacheConfig(
            fitting_results=args.fitting_results,
            annotations=list(args.annotation),
            output=args.output,
            min_size=args.min_size,
            max_size=args.max_size,
            max_term_genes=args.max_term_genes,
            dr_threshold=args.dr_threshold,
            n_permutations=args.n_permutations,
            random_state=args.random_state,
            jobs=args.jobs,
        )
        run(config)
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
