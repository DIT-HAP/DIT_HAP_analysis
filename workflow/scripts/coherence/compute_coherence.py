#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Macromolecular Complex Coherence Analysis
=========================================

Per-dataset: for every macromolecular complex whose DR>0.3 members number
between --min-size and --max-size, measures how tightly its member genes
cluster in the 2D DIT-HAP fitness space and tests that tightness against a
genome-wide null via a seeded permutation test. Ported from
DIT_HAP_pipeline/workflow/notebooks/complex_analysis.ipynb (section 5).

Fitness "points" are the min-max normalized (DR, DL/10) coordinates of each
gene (the notebook's `normalized_um_DITHAP`, `normalized_lam_DITHAP`: DR is
normalized against [0, 1] so it is unchanged, DL against [0, 10] so it is
divided by 10). Coherence = small median pairwise distance (MPD) among members
relative to random draws of the same number of background genes.

Input
-----
- final_clusters.tsv (Systematic ID, A, DR, DL, cluster) from the clustering
  finalize-variant stage; the rule sources it via final_clusters_path(dataset,
  selected_variant). Only Systematic ID / A / DR / DL are read here. Legacy
  releases may still ship the pre-rename um/lam headers -> normalized to DR/DL.
- PomBase macromolecular_complex_annotation.tsv (one row per complex-member:
  complex_term_id, GO_term_name, systematic_id, symbol, ...). Maps complexes ->
  member genes.

Output
------
- complex_coherence_metrics.tsv: one row per surviving complex (complex name,
  term_size, geometric-median centroid, pairwise-distance stats, observed_mpd,
  z_score, p_value, n_permutations).
- coherence_analysis.pdf: complex-size histogram + a coherence volcano
  (z-score vs -log10 permutation p-value, sized by complex size).

Usage
-----
    python compute_complex_coherence.py \\
        --final-clusters results/clustering/final/{dataset}/{variant}/final_clusters.tsv \\
        --complex-annotation .../macromolecular_complex_annotation.tsv \\
        --min-size 3 --max-size 300 --dr-threshold 0.3 \\
        --n-permutations 1000 --random-state 42 \\
        --output-metrics results/complex/{dataset}/complex_coherence_metrics.tsv \\
        --output-figure results/complex/{dataset}/coherence_analysis.pdf

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-20
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
# workflow/src/clustering/candidates.py's _LEGACY_METRIC_RENAME): some curated
# final_clusters.tsv releases still ship the pre-rename um/lam headers.
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}

# Min-max normalization ranges for the DIT-HAP fitness "points", byte-faithful
# to the source notebook: DR against [0, 1] (unchanged), DL against [0, 10]
# (i.e. DL / 10). Points below/above the range are NOT clipped, matching the
# notebook's plain (value - min) / (max - min).
_DR_NORM_RANGE = (0.0, 1.0)
_DL_NORM_RANGE = (0.0, 10.0)

# PomBase annotation column -> canonical name used throughout this script.
_ANNOTATION_RENAME = {"systematic_id": "Systematic ID", "symbol": "Name"}


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class CoherenceConfig:
    """Inputs, outputs, and parameters for the complex coherence analysis."""
    final_clusters: Path
    complex_annotation: Path
    output_metrics: Path
    output_figure: Path
    min_size: int = 3
    max_size: int = 300
    dr_threshold: float = 0.3
    n_permutations: int = 1000
    random_state: int = 42

    def validate(self) -> None:
        """Raise ValueError if inputs are missing or params invalid, then make output dirs."""
        for path in [self.final_clusters, self.complex_annotation]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.min_size < 2:
            raise ValueError(f"min_size must be >= 2 (pairwise distances need 2 points): {self.min_size}")
        if self.max_size < self.min_size:
            raise ValueError(f"max_size ({self.max_size}) must be >= min_size ({self.min_size})")
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


def load_final_clusters(final_clusters_path: Path, dr_threshold: float) -> pd.DataFrame:
    """Load curated final_clusters.tsv, normalize legacy um/lam -> DR/DL, add fitness points.

    Adds `norm_DR`/`norm_DL` (the 2D coherence coordinates) and keeps only the
    DR>dr_threshold, fully-fitted genes that make up the genome-wide background.
    """
    clusters = pd.read_csv(final_clusters_path, sep="\t")
    rename = {
        old: new
        for old, new in _LEGACY_METRIC_RENAME.items()
        if old in clusters.columns and new not in clusters.columns
    }
    if rename:
        logger.info(f"Normalizing legacy metric columns: {rename}")
        clusters = clusters.rename(columns=rename)

    for required in ["Systematic ID", "DR", "DL"]:
        if required not in clusters.columns:
            raise ValueError(f"final_clusters.tsv missing required column '{required}' (have: {list(clusters.columns)})")

    # Genes must have both FINITE fitness coordinates to be placed in the 2D
    # space. Map +/-inf to NaN first so dropna removes it too: an inf DR/DL
    # would normalize to inf and poison pdist (observed_mpd=nan, z_score=nan,
    # and `null_mpds <= nan` -> all-False -> a spurious p_value=0.0 row).
    clusters = clusters.replace([np.inf, -np.inf], np.nan).dropna(subset=["DR", "DL"]).copy()
    clusters["norm_DR"] = _min_max_normalize(clusters["DR"].to_numpy(dtype=float), *_DR_NORM_RANGE)
    clusters["norm_DL"] = _min_max_normalize(clusters["DL"].to_numpy(dtype=float), *_DL_NORM_RANGE)

    # Section 5.2: filter to the non-WT / depleting genes (DR > threshold).
    background = clusters[clusters["DR"] > dr_threshold].copy()
    logger.info(
        f"final_clusters.tsv: {len(clusters):,} fitted genes -> "
        f"{len(background):,} background genes with DR > {dr_threshold}"
    )
    return background


def load_complex_annotation(annotation_path: Path) -> pd.DataFrame:
    """Load PomBase macromolecular_complex_annotation.tsv (complex -> member genes)."""
    annotation = pd.read_csv(annotation_path, sep="\t").rename(columns=_ANNOTATION_RENAME)
    for required in ["complex_term_id", "GO_term_name", "Systematic ID"]:
        if required not in annotation.columns:
            raise ValueError(
                f"complex annotation missing required column '{required}' (have: {list(annotation.columns)})"
            )
    keep = [c for c in ["complex_term_id", "GO_term_name", "Systematic ID", "Name"] if c in annotation.columns]
    return annotation[keep].drop_duplicates()


# =============================================================================
# CORE LOGIC — coherence per complex
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


def build_complex_groups(
    background: pd.DataFrame, annotation: pd.DataFrame, min_size: int, max_size: int
) -> dict[str, pd.DataFrame]:
    """Map surviving complexes -> their DR>threshold member rows.

    Merges the DR>threshold background genes onto the complex annotation
    (inner join, so only members present in the background survive), groups by
    complex name, and keeps groups whose member count is within
    [min_size, max_size]. Byte-faithful to the notebook's section 5.2/5.3:
    the size filter counts DR>threshold members, and the notebook itself keys
    coherence groups on `GO_term_name` (the human-readable complex name).
    """
    # Merge complex membership onto the already DR-filtered background so both
    # the observed complex points AND the eligible member set are DR>threshold.
    merged = annotation.merge(
        background[["Systematic ID", "norm_DR", "norm_DL"]], on="Systematic ID", how="inner"
    )
    groups = {}
    for name, group in merged.groupby("GO_term_name"):
        # A gene may be annotated to a complex more than once; dedupe by gene.
        group = group.drop_duplicates(subset="Systematic ID")
        if min_size <= len(group) <= max_size:
            groups[name] = group
    logger.info(
        f"{merged['GO_term_name'].nunique():,} complexes with >=1 background member -> "
        f"{len(groups):,} with {min_size} <= size <= {max_size}"
    )
    return groups


def compute_coherence_table(
    groups: dict[str, pd.DataFrame],
    background_points: np.ndarray,
    background_index: dict[str, int],
    n_permutations: int,
    random_state: int,
) -> pd.DataFrame:
    """One coherence row per complex: metrics + permutation z-score of the MPD.

    `background_points` is the (n_background, 2) genome-wide point cloud;
    `background_index` maps Systematic ID -> row index into it, so a complex's
    members are addressed as row indices for the permutation null (the null
    draws random background rows of the same count).
    """
    rows = []
    for name, group in groups.items():
        member_ids = group["Systematic ID"].tolist()
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
            "complex": name,
            "complex_term_id": group["complex_term_id"].iloc[0],
            "term_size": len(group),
            "covered_genes": ", ".join(sorted(group["Name"].dropna().astype(str))) if "Name" in group else "",
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
    """Load -> filter -> per-complex coherence + permutation test -> TSV + figure."""
    config.validate()

    background = load_final_clusters(config.final_clusters, config.dr_threshold)
    annotation = load_complex_annotation(config.complex_annotation)

    # Genome-wide background point cloud + Systematic ID -> row index map.
    background = background.reset_index(drop=True)
    background_points = background[["norm_DR", "norm_DL"]].to_numpy(dtype=float)
    background_index = {gid: i for i, gid in enumerate(background["Systematic ID"])}

    groups = build_complex_groups(background, annotation, config.min_size, config.max_size)
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
        f"Coherence: {len(table):,} complexes scored, {n_coherent:,} coherent (z<0); "
        f"wrote {config.output_metrics}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Compute macromolecular complex coherence metrics")
    parser.add_argument("--final-clusters", type=Path, required=True, help="Curated final_clusters.tsv")
    parser.add_argument("--complex-annotation", type=Path, required=True, help="PomBase macromolecular_complex_annotation.tsv")
    parser.add_argument("--min-size", type=int, default=3, help="Minimum DR>threshold members per complex")
    parser.add_argument("--max-size", type=int, default=300, help="Maximum DR>threshold members per complex")
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
            final_clusters=args.final_clusters,
            complex_annotation=args.complex_annotation,
            output_metrics=args.output_metrics,
            output_figure=args.output_figure,
            min_size=args.min_size,
            max_size=args.max_size,
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
