#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gene-Group Coherence Analysis (source-agnostic)
===============================================

Per-dataset x source: for every group (complex / GO term / ...) whose
DR>threshold members number between --min-size and --max-size AND whose total
annotated membership is <= --max-term-genes, measures how tightly its member
genes cluster in the 2D DIT-HAP fitness space and tests that tightness against
a genome-wide null via a seeded permutation test. Ported from
DIT_HAP_pipeline/workflow/notebooks/complex_analysis.ipynb (section 5); now
generalized to any grouping source via a prepared long-table.

Fitness "points" are the min-max normalized (DR, DL/10) coordinates of each
gene (the notebook's `normalized_um_DITHAP`, `normalized_lam_DITHAP`: DR is
normalized against [0, 1] so it is unchanged, DL against [0, 10] so it is
divided by 10). Coherence = small median pairwise distance (MPD) among members
relative to random draws of the same number of background genes.

Input
-----
- fitting_results.tsv: the upstream per-gene fitting statistics, systematic id
  as the index (column 0), with DR/DL fitness columns. Legacy releases may still
  ship the pre-rename um/lam headers -> normalized to DR/DL. Only the systematic
  id + DR/DL are read here.
- group_annotation_long.tsv: the prepared unified long-table from
  prepare_annotation.py (one row per group-member), with the contract columns
  (source, group_id, group_name, Systematic ID, Name, n_group_genes). Maps
  groups -> member genes. n_group_genes = total annotated members (pre-DR-filter).

Output
------
- coherence_metrics.tsv: one row per surviving group, in emission order:
  source, group_id, group_name, term_size, n_group_genes, covered_genes, then
  the geometric-median centroid + pairwise-distance stats (centroid_x,
  centroid_y, median/mean/std/min/max_distance, mpd), observed_mpd, z_score,
  p_value, n_permutations.
- coherence_analysis.pdf: group-size histogram + a coherence volcano
  (z-score vs -log10 permutation p-value, sized by group size).

Usage
-----
    python compute_coherence.py \\
        --fitting-results .../fitting_results.tsv \\
        --annotation results/coherence/{dataset}/{source}/group_annotation_long.tsv \\
        --source go_macrocomplex \\
        --min-size 3 --max-size 300 --max-term-genes 500 --dr-threshold 0.3 \\
        --n-permutations 1000 --random-state 42 \\
        --output-metrics results/coherence/{dataset}/{source}/coherence_metrics.tsv \\
        --output-figure results/coherence/{dataset}/{source}/coherence_analysis.pdf

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-20
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

# 2. Data Processing Imports
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

# 3. Third-party Imports
import matplotlib

matplotlib.use("Agg")  # headless: this script only writes a PDF, never displays
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
# The Weiszfeld geometric median + the seeded MPD permutation test are the
# shared coherence algorithm, sourced from the canonical workflow/src/coherence
# module (also used by the themes-A/D verification_complex scripts). We call it
# with median_pairwise_distance, matching this analysis's coherence axis; the
# descriptive pairwise-distance summary below is trivial numpy kept local (the
# shared coherence_metrics exposes a different, richer key set).
from workflow.src.coherence.metrics import (  # noqa: E402
    geometric_median,
    compute_distance_zscore,
)
from workflow.src.plotting.style import AX_HEIGHT, AX_WIDTH  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Legacy -> current metric column names (same quirk as
# workflow/src/clustering/candidates.py's _LEGACY_METRIC_RENAME): some upstream
# fitting_results.tsv exports still ship the pre-rename um/lam headers.
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}

# Min-max normalization ranges for the DIT-HAP fitness "points", byte-faithful
# to the source notebook: DR against [0, 1] (unchanged), DL against [0, 10]
# (i.e. DL / 10). Points below/above the range are NOT clipped, matching the
# notebook's plain (value - min) / (max - min).
_DR_NORM_RANGE = (0.0, 1.0)
_DL_NORM_RANGE = (0.0, 10.0)

# The prepared long-table contract (from prepare_annotation.py / sources.py).
_LONG_TABLE_COLUMNS = [
    "source", "group_id", "group_name", "Systematic ID", "Name", "n_group_genes",
]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class CoherenceConfig:
    """Inputs, outputs, and parameters for the gene-group coherence analysis."""
    fitting_results: Path
    annotation: Path
    source: str
    output_metrics: Path
    output_figure: Path
    min_size: int = 3
    max_size: int = 300
    max_term_genes: int = 500
    dr_threshold: float = 0.3
    n_permutations: int = 1000
    random_state: int = 42
    features: Path | None = None  # WIRED but unused here; Task 5's figure consumes it.

    def validate(self) -> None:
        """Raise ValueError if inputs are missing or params invalid, then make output dirs."""
        required = [self.fitting_results, self.annotation]
        if self.features is not None:
            required.append(self.features)
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
        for out in [self.output_metrics, self.output_figure]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


def _min_max_normalize(values: np.ndarray, min_value: float, max_value: float) -> np.ndarray:
    """Normalize to [0, 1] via (v - min) / (max - min); all-zeros if the range is degenerate.

    Byte-faithful to the notebook's min_max_normalization: values outside
    [min_value, max_value] are NOT clipped, so a DR just above 1.0 maps just
    above 1.0. Only the max==min degenerate case is guarded.
    """
    if max_value - min_value == 0:
        return np.zeros_like(values, dtype=float)
    return (values - min_value) / (max_value - min_value)


def load_fitting_results(fitting_results_path: Path, dr_threshold: float) -> pd.DataFrame:
    """Load upstream fitting_results.tsv, normalize legacy um/lam -> DR/DL, add fitness points.

    The systematic id is the INDEX (column 0), matching how
    workflow/src/clustering/candidates.py reads it; we reset it to a
    `Systematic ID` column (robust to whatever the index was named). Adds
    `norm_DR`/`norm_DL` (the 2D coherence coordinates) and keeps only the
    DR>dr_threshold, fully-fitted genes that make up the genome-wide background.
    """
    fitting = pd.read_csv(fitting_results_path, sep="\t", index_col=0)
    fitting = fitting.reset_index()
    # After reset_index the id column keeps the index's name (e.g. "Systematic ID",
    # "gene_systematic_id", or "index" if the index was unnamed). Normalize it.
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

    # Genes must have both FINITE fitness coordinates to be placed in the 2D
    # space. Map +/-inf to NaN first so dropna removes it too: an inf DR/DL
    # would normalize to inf and poison pdist (observed_mpd=nan, z_score=nan,
    # and `null_mpds <= nan` -> all-False -> a spurious p_value=0.0 row).
    fitting = fitting.replace([np.inf, -np.inf], np.nan).dropna(subset=["DR", "DL"]).copy()
    fitting["norm_DR"] = _min_max_normalize(fitting["DR"].to_numpy(dtype=float), *_DR_NORM_RANGE)
    fitting["norm_DL"] = _min_max_normalize(fitting["DL"].to_numpy(dtype=float), *_DL_NORM_RANGE)

    # Section 5.2: filter to the non-WT / depleting genes (DR > threshold).
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
    """Geometric-median centroid + descriptive stats of all pairwise L2 distances.

    The centroid uses the shared Weiszfeld geometric_median; the distance
    summary is a trivial numpy reduction over pdist kept local so this script's
    output-TSV schema (median/mean/std/min/max_distance + mpd) is stable and
    independent of the shared coherence_metrics' richer key set. `mpd` (median
    pairwise distance) is the coherence statistic the permutation test scores
    and equals `median_distance`. Degenerate inputs (0 or 1 point) have no
    pairwise distances; their stats are reported as 0.0.
    """
    points = np.asarray(points, dtype=float)
    centroid = geometric_median(points)

    pairwise = pdist(points)  # condensed vector of all C(n,2) L2 distances
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
        "mpd": median_d,  # median pairwise distance == median_distance
    }


def build_groups(
    background: pd.DataFrame,
    long_table: pd.DataFrame,
    min_group_size: int,
    max_group_size: int,
    max_term_genes: int,
) -> dict[str, pd.DataFrame]:
    """Map surviving groups (keyed on group_id) -> their DR>threshold member rows.

    Merges the DR>threshold background genes onto the long-table (inner join, so
    only members present in the background survive), groups by `group_id`, and
    keeps groups that pass BOTH filters:
      - total annotated membership (`n_group_genes`, pre-DR-filter) <= max_term_genes,
        so overly-broad terms are dropped before scoring; and
      - DR>threshold member count within [min_group_size, max_group_size].
    Byte-faithful to the notebook's section 5.2/5.3 (the size filter counts
    DR>threshold members); keying on the stable `group_id` (not the human name)
    matches how the long-table dedups/counts in sources.py.
    """
    # Merge group membership onto the already DR-filtered background so both the
    # observed group points AND the eligible member set are DR>threshold.
    merged = long_table.merge(
        background[["Systematic ID", "norm_DR", "norm_DL"]], on="Systematic ID", how="inner"
    )
    groups = {}
    for group_id, grp in merged.groupby("group_id"):
        # A gene may be annotated to a group more than once; dedupe by gene.
        grp = grp.drop_duplicates(subset="Systematic ID")
        # n_group_genes is per-group-constant (sources.py sets it via
        # groupby(group_id).transform("size")), so .iloc[0] safely reads it.
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
    """One coherence row per group: identity + metrics + permutation z-score of the MPD.

    `background_points` is the (n_background, 2) genome-wide point cloud;
    `background_index` maps Systematic ID -> row index into it, so a group's
    members are addressed as row indices for the permutation null (the null
    draws random background rows of the same count).
    """
    rows = []
    for group_id, grp in groups.items():
        member_ids = grp["Systematic ID"].tolist()
        member_indices = [background_index[gid] for gid in member_ids]
        member_points = background_points[member_indices]

        metrics = coherence_metrics(member_points)
        # Shared permutation test: X = member points, bg = the FULL background
        # point cloud (members included, matching the notebook's null draw), and
        # median_pairwise_distance is this analysis's coherence axis. Returns a
        # (z_score, p_value) tuple; observed_mpd is the local metrics["mpd"].
        z_score, p_value = compute_distance_zscore(
            member_points,
            background_points,
            method="median_pairwise_distance",
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
            "observed_mpd": metrics["mpd"],
            "z_score": z_score,
            "p_value": p_value,
            "n_permutations": n_permutations,
        })

    table = pd.DataFrame(rows)
    if not table.empty:
        table = table.sort_values("z_score").reset_index(drop=True)
    return table


# =============================================================================
# PLOTTING
# =============================================================================
def plot_coherence(table: pd.DataFrame) -> plt.Figure:
    """Complex-size histogram + coherence volcano (z-score vs -log10 p-value)."""
    fig, axes = plt.subplots(1, 2, figsize=(AX_WIDTH * 2, AX_HEIGHT))

    ax_size, ax_volcano = axes
    if table.empty:
        for ax in axes:
            ax.text(0.5, 0.5, "No complexes passed the size filter", ha="center", va="center")
        fig.tight_layout()
        return fig

    ax_size.hist(table["term_size"], bins=20, rwidth=0.9, color="#6b99df")
    ax_size.set_xlabel("Complex size (DR>threshold members)")
    ax_size.set_ylabel("Number of complexes")
    ax_size.set_title("Complex size distribution")

    # Volcano: negative z-score (tight/coherent) to the left, significance up.
    # Guard p==0 (observed below every permutation) so -log10 stays finite by
    # flooring at the smallest resolvable p-value, 1 / n_permutations.
    n_perm = int(table["n_permutations"].iloc[0])
    p_floor = 1.0 / n_perm
    neglog_p = -np.log10(table["p_value"].clip(lower=p_floor))
    sizes = np.sqrt(table["term_size"].to_numpy(dtype=float)) * 8
    scatter = ax_volcano.scatter(
        table["z_score"], neglog_p, s=sizes, c=table["z_score"], cmap="coolwarm_r", alpha=0.8, edgecolors="none"
    )
    ax_volcano.axvline(0.0, color="gray", linestyle="--", linewidth=0.8)
    ax_volcano.set_xlabel("MPD z-score (negative = coherent)")
    ax_volcano.set_ylabel("-log10(permutation p-value)")
    ax_volcano.set_title("Complex coherence")
    fig.colorbar(scatter, ax=ax_volcano, label="z-score")

    fig.tight_layout()
    return fig


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: CoherenceConfig) -> None:
    """Load -> filter -> per-group coherence + permutation test -> TSV + figure."""
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
    table.to_csv(config.output_metrics, sep="\t", index=False)

    fig = plot_coherence(table)
    with PdfPages(config.output_figure) as pdf:
        pdf.savefig(fig, dpi=300, bbox_inches="tight")
    plt.close(fig)

    n_coherent = int((table["z_score"] < 0).sum()) if not table.empty else 0
    logger.success(
        f"[{config.source}] Coherence: {len(table):,} groups scored, {n_coherent:,} coherent (z<0); "
        f"wrote {config.output_metrics}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Compute gene-group coherence metrics")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Upstream fitting_results.tsv (systematic id as index col 0)")
    parser.add_argument("--annotation", type=Path, required=True, help="Prepared group_annotation_long.tsv (long-table from prepare_annotation.py)")
    parser.add_argument("--source", type=str, required=True, help="Grouping-database source name (fan-out dimension)")
    parser.add_argument("--features", type=Path, default=None, help="Optional gene features TSV (wired for Task 5's figure; unused here)")
    parser.add_argument("--min-size", type=int, default=3, help="Minimum DR>threshold members per group")
    parser.add_argument("--max-size", type=int, default=300, help="Maximum DR>threshold members per group")
    parser.add_argument("--max-term-genes", type=int, default=500, help="Drop groups whose total annotated membership (n_group_genes) exceeds this")
    parser.add_argument("--dr-threshold", type=float, default=0.3, help="Keep genes with DR > this")
    parser.add_argument("--n-permutations", type=int, default=1000, help="Permutation null draws")
    parser.add_argument("--random-state", type=int, default=42, help="Permutation RNG seed")
    parser.add_argument("--output-metrics", type=Path, required=True, help="Output coherence metrics TSV")
    parser.add_argument("--output-figure", type=Path, required=True, help="Output coherence figure PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
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
            output_metrics=args.output_metrics,
            output_figure=args.output_figure,
            min_size=args.min_size,
            max_size=args.max_size,
            max_term_genes=args.max_term_genes,
            dr_threshold=args.dr_threshold,
            n_permutations=args.n_permutations,
            random_state=args.random_state,
            features=args.features,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
