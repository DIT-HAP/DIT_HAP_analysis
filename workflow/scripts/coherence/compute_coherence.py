#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gene-Group Coherence Analysis (source-agnostic) — Computation Only
==================================================================

Per-dataset x source: for every group (complex / GO term / ...) whose
DR>threshold members number between --min-size and --max-size AND whose total
annotated membership is <= --max-term-genes, measures how tightly its member
genes cluster in the 2D DIT-HAP fitness space and tests that tightness against
a genome-wide null via a seeded permutation test.

This is the computation component after the computation-plotting decoupling
refactor (ADR-0001). The visualization is handled by plot_coherence.py.

Fitness "points" are the min-max normalized (DR, DL/10) coordinates of each
gene. Coherence = small median pairwise distance (MPD) among members relative
to random draws of the same number of background genes.

Input
-----
- fitting_results.tsv: the upstream per-gene fitting statistics, systematic id
  as the index (column 0), with DR/DL fitness columns. Legacy releases may still
  ship the pre-rename um/lam headers -> normalized to DR/DL.
- group_annotation_long.tsv: the prepared unified long-table from
  prepare_annotation.py (one row per group-member), with the contract columns
  (source, group_id, group_name, Systematic ID, Name, n_group_genes).

Output
------
- coherence.parquet: one row per surviving group with columns: source, group_id,
  group_name, term_size, n_group_genes, covered_genes, centroid_x, centroid_y,
  median_distance, mean_distance, std_distance, min_distance, max_distance, mpd,
  z_score, p_value, n_permutations, mean_pairwise_distance_zscore,
  mean_pairwise_distance_p_value, p_fdr.

Usage
-----
    python compute_coherence.py \\
        --fitting-results .../fitting_results.tsv \\
        --annotation results/coherence/{dataset}/{source}/group_annotation_long.tsv \\
        --source go_macrocomplex \\
        --min-size 3 --max-size 300 --max-term-genes 500 --dr-threshold 0.3 \\
        --n-permutations 1000 --random-state 42 \\
        --output results/coherence/{dataset}/{source}/coherence.parquet

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (refactor)
Date:     2026-09-03
Version:  3.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Third-party
import numpy as np
import pandas as pd
from loguru import logger
from scipy.spatial.distance import pdist
from scipy.stats import false_discovery_control

# 3. Local
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from coherence.metrics import (  # noqa: E402
    geometric_median,
    compute_distance_zscore,
)
from io_table import write_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Legacy -> current metric column names
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}

# Min-max normalization ranges for the DIT-HAP fitness "points"
_DR_NORM_RANGE = (0.0, 1.0)
_DL_NORM_RANGE = (0.0, 10.0)

# The prepared long-table contract
_LONG_TABLE_COLUMNS = [
    "source", "group_id", "group_name", "Systematic ID", "Name", "n_group_genes",
]


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class CoherenceConfig:
    """Inputs, outputs, and parameters for the gene-group coherence analysis."""
    fitting_results: Path
    annotation: Path
    source: str
    output: Path
    min_size: int = 3
    max_size: int = 300
    max_term_genes: int = 500
    dr_threshold: float = 0.3
    n_permutations: int = 1000
    random_state: int = 42

    def validate(self) -> None:
        """Raise ValueError if inputs are missing or params invalid, then make output dirs."""
        required = [self.fitting_results, self.annotation]
        for path in required:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.min_size < 2:
            raise ValueError(f"min_size must be >= 2 (pairwise distances need 2 points): {self.min_size}")
        if self.max_size < self.min_size:
            raise ValueError(f"max_size ({self.max_size}) must be >= min_size ({self.min_size})")
        if self.max_term_genes < self.min_size:
            raise ValueError(f"max_term_genes ({self.max_term_genes}) must be >= min_size ({self.min_size})")
        if self.n_permutations < 1:
            raise ValueError(f"n_permutations must be >= 1: {self.n_permutations}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def _min_max_normalize(values: np.ndarray, min_value: float, max_value: float) -> np.ndarray:
    """Normalize to [0, 1] via (v - min) / (max - min); all-zeros if the range is degenerate."""
    if max_value - min_value == 0:
        return np.zeros_like(values, dtype=float)
    return (values - min_value) / (max_value - min_value)


def load_fitting_results(fitting_results_path: Path, dr_threshold: float) -> pd.DataFrame:
    """Load upstream fitting_results.tsv, normalize legacy um/lam -> DR/DL, add fitness points."""
    fitting = pd.read_csv(fitting_results_path, sep="\t", index_col=0)
    fitting = fitting.reset_index()

    if "Systematic ID" not in fitting.columns:
        first_col = fitting.columns[0]
        logger.info(f"Renaming fitting_results index column '{first_col}' -> 'Systematic ID'")
        fitting = fitting.rename(columns={first_col: "Systematic ID"})

    rename = {
        old: new
        for old, new in _LEGACY_METRIC_RENAME.items()
        if old in fitting.columns and new not in fitting.columns
    }
    if rename:
        logger.info(f"Normalizing legacy metric columns: {rename}")
        fitting = fitting.rename(columns=rename)

    for required in ["Systematic ID", "DR", "DL"]:
        if required not in fitting.columns:
            raise ValueError(f"fitting_results.tsv missing required column '{required}' (have: {list(fitting.columns)})")

    fitting = fitting.replace([np.inf, -np.inf], np.nan).dropna(subset=["DR", "DL"]).copy()
    fitting["norm_DR"] = _min_max_normalize(fitting["DR"].to_numpy(dtype=float), *_DR_NORM_RANGE)
    fitting["norm_DL"] = _min_max_normalize(fitting["DL"].to_numpy(dtype=float), *_DL_NORM_RANGE)

    background = fitting[fitting["DR"] > dr_threshold].copy()
    logger.info(
        f"fitting_results.tsv: {len(fitting):,} fitted genes -> "
        f"{len(background):,} background genes with DR > {dr_threshold}"
    )
    return background


def load_long_table(annotation_path: Path) -> pd.DataFrame:
    """Load the prepared unified long-table (group -> member genes) + validate its contract."""
    long_table = pd.read_csv(annotation_path, sep="\t")
    for required in _LONG_TABLE_COLUMNS:
        if required not in long_table.columns:
            raise ValueError(
                f"annotation long-table missing required column '{required}' (have: {list(long_table.columns)})"
            )
    return long_table


# =============================================================================
# CORE LOGIC — coherence per group
# =============================================================================
def coherence_metrics(points: np.ndarray) -> dict:
    """Geometric-median centroid + descriptive stats of all pairwise L2 distances."""
    points = np.asarray(points, dtype=float)
    centroid = geometric_median(points)

    pairwise = pdist(points)
    if pairwise.size == 0:
        median_d = mean_d = std_d = min_d = max_d = 0.0
    else:
        median_d = float(np.median(pairwise))
        mean_d = float(np.mean(pairwise))
        std_d = float(np.std(pairwise))
        min_d = float(np.min(pairwise))
        max_d = float(np.max(pairwise))

    return {
        "centroid_x": float(centroid[0]),
        "centroid_y": float(centroid[1]),
        "median_distance": median_d,
        "mean_distance": mean_d,
        "std_distance": std_d,
        "min_distance": min_d,
        "max_distance": max_d,
        "mpd": median_d,
    }


def build_groups(
    background: pd.DataFrame,
    long_table: pd.DataFrame,
    min_group_size: int,
    max_group_size: int,
    max_term_genes: int,
) -> dict[str, pd.DataFrame]:
    """Map surviving groups (keyed on group_id) -> their DR>threshold member rows."""
    merged = long_table.merge(
        background[["Systematic ID", "norm_DR", "norm_DL"]], on="Systematic ID", how="inner"
    )
    groups = {}
    for group_id, grp in merged.groupby("group_id"):
        grp = grp.drop_duplicates(subset="Systematic ID")
        n_total = int(grp["n_group_genes"].iloc[0])
        if n_total > max_term_genes:
            continue
        if min_group_size <= len(grp) <= max_group_size:
            groups[group_id] = grp
    logger.info(
        f"{merged['group_id'].nunique():,} groups with >=1 background member -> "
        f"{len(groups):,} with {min_group_size} <= size <= {max_group_size} "
        f"and n_group_genes <= {max_term_genes}"
    )
    return groups


def compute_coherence_table(
    groups: dict[str, pd.DataFrame],
    background_points: np.ndarray,
    background_index: dict[str, int],
    n_permutations: int,
    random_state: int,
) -> pd.DataFrame:
    """One coherence row per group: identity + metrics + permutation z-score of the MPD."""
    rows = []
    for group_id, grp in groups.items():
        member_ids = grp["Systematic ID"].tolist()
        member_indices = [background_index[gid] for gid in member_ids]
        member_points = background_points[member_indices]

        metrics = coherence_metrics(member_points)
        z_score, p_value = compute_distance_zscore(
            member_points,
            background_points,
            method="median_pairwise_distance",
            n_permutations=n_permutations,
            random_state=random_state,
        )

        mean_pairwise_distance_zscore, mean_pairwise_distance_p_value = compute_distance_zscore(
            member_points,
            background_points,
            method="mean_pairwise_distance",
            n_permutations=n_permutations,
            random_state=random_state,
        )

        rows.append({
            "source": grp["source"].iloc[0],
            "group_id": group_id,
            "group_name": grp["group_name"].iloc[0],
            "term_size": len(grp),
            "n_group_genes": int(grp["n_group_genes"].iloc[0]),
            "covered_genes": ", ".join(sorted(grp["Name"].dropna().astype(str))) if "Name" in grp else "",
            **metrics,
            "z_score": z_score,
            "p_value": p_value,
            "n_permutations": n_permutations,
            "mean_pairwise_distance_zscore": mean_pairwise_distance_zscore,
            "mean_pairwise_distance_p_value": mean_pairwise_distance_p_value,
        })

    table = pd.DataFrame(rows)
    if not table.empty:
        table["p_fdr"] = false_discovery_control(table["p_value"].to_numpy(), method="bh")
        table = table.sort_values("z_score").reset_index(drop=True)
    return table


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: CoherenceConfig) -> None:
    """Load -> filter -> per-group coherence + permutation test -> Parquet."""
    config.validate()

    background = load_fitting_results(config.fitting_results, config.dr_threshold)
    long_table = load_long_table(config.annotation)

    # Genome-wide background point cloud + Systematic ID -> row index map.
    background = background.reset_index(drop=True)
    background_points = background[["norm_DR", "norm_DL"]].to_numpy(dtype=float)
    background_index = {gid: i for i, gid in enumerate(background["Systematic ID"])}

    groups = build_groups(
        background, long_table, config.min_size, config.max_size, config.max_term_genes
    )
    table = compute_coherence_table(
        groups, background_points, background_index, config.n_permutations, config.random_state
    )

    # Write output as Parquet
    write_parquet(table, config.output)

    n_coherent = int((table["z_score"] < 0).sum()) if not table.empty else 0
    logger.success(
        f"[{config.source}] Coherence: {len(table):,} groups scored, {n_coherent:,} coherent (z<0); "
        f"wrote {config.output}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Compute gene-group coherence metrics")
    parser.add_argument("--fitting-results", type=Path, required=True,
                        help="Upstream fitting_results.tsv (systematic id as index col 0)")
    parser.add_argument("--annotation", type=Path, required=True,
                        help="Prepared group_annotation_long.tsv (long-table from prepare_annotation.py)")
    parser.add_argument("--source", type=str, required=True,
                        help="Grouping-database source name (fan-out dimension)")
    parser.add_argument("--min-size", type=int, default=3,
                        help="Minimum DR>threshold members per group")
    parser.add_argument("--max-size", type=int, default=300,
                        help="Maximum DR>threshold members per group")
    parser.add_argument("--max-term-genes", type=int, default=500,
                        help="Drop groups whose total annotated membership (n_group_genes) exceeds this")
    parser.add_argument("--dr-threshold", type=float, default=0.3,
                        help="Keep genes with DR > this")
    parser.add_argument("--n-permutations", type=int, default=1000,
                        help="Permutation null draws")
    parser.add_argument("--random-state", type=int, default=42,
                        help="Permutation RNG seed")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output coherence metrics Parquet")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = CoherenceConfig(
            fitting_results=args.fitting_results,
            annotation=args.annotation,
            source=args.source,
            output=args.output,
            min_size=args.min_size,
            max_size=args.max_size,
            max_term_genes=args.max_term_genes,
            dr_threshold=args.dr_threshold,
            n_permutations=args.n_permutations,
            random_state=args.random_state,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
