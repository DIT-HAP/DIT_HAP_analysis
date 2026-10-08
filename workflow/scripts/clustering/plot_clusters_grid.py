#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Clusters Grid — one PDF per k, all methods side by side
============================================================

Renders the (method x k) label table from compute_clustering_grid as one PDF
per k: page is a grid of DR/DL scatters, one subplot per method, all on shared
axes so a reader can see how the methods partition the same space at the same
cluster count (the visual companion to k_selection.pdf's numeric QC). Cluster
colors are one shared mapping across every page (DR-ascending ids), so the
same cluster reads as the same colour everywhere.

At the reference k (config grid_finalize.k) an extra panel renders the
HAND-CURATED cluster annotation (manual_curated_clusters.tsv) under the id
"manual", on the same axes — the human judgment the numeric metrics cannot
score, placed next to the methods it should be compared against.

Input
-----
- grid_labels.parquet (from compute_clustering_grid)
- scaled_data.parquet (the DR/DL coordinates to scatter)
- manual_curated_clusters.tsv (only for the reference-k panel)

Output
------
- clusters_grid_k<k>.pdf + .review.png sibling, one per k in the table.

Usage
-----
    python plot_clusters_grid.py \\
        --labels results/clustering/{dataset}/_work/grid_labels.parquet \\
        --scaled results/clustering/{dataset}/_work/scaled_data.parquet \\
        --output-stem results/clustering/{dataset}/clusters_grid \\
        --dataset HD_DIT_HAP --reference-k 9 \\
        --curated tmp/manual_curated_clusters.tsv

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-10-08
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
import matplotlib

matplotlib.use("Agg")  # headless: writes PDF/PNG only
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from figures import apply_house_style, grid_axes, save_dual  # noqa: E402
from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402
from plotting.gene_level import CLUSTER_COLORS  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotGridConfig:
    """Label table + coordinates in; one PDF per k out (stem + _k<k>)."""
    labels: Path
    scaled: Path
    annotated: Path
    output_stem: Path
    dataset: str
    reference_k: int | None = None
    curated: Path | None = None
    col_num: int = 3

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [self.labels, self.scaled, self.annotated]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.reference_k is not None and not self.curated.exists():
            raise ValueError(f"--curated table not found: {self.curated} (required for the reference-k panel)")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def scatter_clusters(ax: plt.Axes, df: pd.DataFrame, colors: dict[int, str]) -> None:
    """One method's DR/DL scatter on one axes: per-cluster color + centroid id chip."""
    for cluster, cluster_df in df.groupby("cluster"):
        ax.scatter(
            cluster_df["DR"], cluster_df["DL"], c=colors[int(cluster)], s=3, linewidths=0, rasterized=True
        )
        ax.annotate(
            f"{cluster}",
            (cluster_df["DR"].mean(), cluster_df["DL"].mean()),
            fontweight="bold",
            fontsize=6,
            ha="center",
            va="center",
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.8),
        )
    ax.set_xlabel("DR")
    ax.set_ylabel("DL")


def curated_labels(curated: Path, annotated: pd.DataFrame) -> pd.Series:
    """The hand-curated clusters, DR-ascending-renumbered to 1..k like every method."""
    from clustering.candidates import renumber_by_dr

    table = pd.read_csv(curated, sep="\t").set_index("Systematic ID")["cluster"]
    labels = table.reindex(annotated.index).dropna().astype(int)
    k = int(labels.nunique())
    return renumber_by_dr(annotated, labels, n_clusters=k, wt_cluster=k)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotGridConfig) -> None:
    """Render one PDF per k: every method side by side (+ the curated panel at the reference k)."""
    config.validate()

    labels = pd.read_parquet(config.labels)
    scaled = pd.read_parquet(config.scaled)
    methods = sorted(labels["method"].unique())

    # One color per DR-ascending cluster id across ALL pages: id c -> house color.
    # Colours beyond the palette wrap (large k): ids stay distinguishable by their
    # centroid chip even when a colour repeats.
    max_cluster = int(labels["cluster"].max())
    colors = {c: CLUSTER_COLORS[(c - 1) % len(CLUSTER_COLORS)] for c in range(1, max_cluster + 1)}

    apply_house_style()
    panel_names = list(methods) + (["manual"] if config.reference_k is not None else [])
    col_num = min(config.col_num, len(panel_names))
    row_num = -(-len(panel_names) // col_num)  # ceiling division

    # The curated renumbering runs on RAW DR/DL (same rule as the methods, same
    # inputs), which is what annotated_data.parquet holds — an extra spine input
    # rather than an un-scale of the clamped matrix (lossy, wrong order risk).
    annotated = read_parquet(config.annotated) if config.reference_k is not None else None
    hand = curated_labels(config.curated, annotated) if config.reference_k is not None else None

    for k in sorted(labels["k"].unique()):
        names = panel_names if k == config.reference_k else list(methods)
        axes = grid_axes(row_num, col_num, labels=None)
        fig = plt.gcf()
        fig.subplots_adjust(top=0.88, hspace=0.45, wspace=0.25)
        for name, ax in zip(names, axes):
            if name == "manual":
                plot_df = scaled.loc[hand.index].assign(cluster=hand.to_numpy())
            else:
                df = labels[(labels["method"] == name) & (labels["k"] == k)]
                plot_df = scaled.loc[df["gene"]].assign(cluster=df["cluster"].to_numpy())
            scatter_clusters(ax, plot_df, colors)
            ax.set_title(f"{name} (k={k}, n={len(plot_df)})", fontsize=7)
        for ax in axes[len(names):]:
            ax.set_visible(False)
        fig.suptitle(f"{config.dataset} — all methods at k={k}", fontsize=10)
        save_dual(config.output_stem.parent / f"{config.output_stem.name}_k{k}")
        plt.close(fig)
    logger.success(f"Wrote {labels['k'].nunique()} per-k PDFs -> {config.output_stem}_k*.pdf")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Plot one PDF per k with all methods side by side")
    parser.add_argument("--labels", type=Path, required=True, help="grid_labels.parquet from compute_clustering_grid")
    parser.add_argument("--scaled", type=Path, required=True, help="scaled_data.parquet from the spine")
    parser.add_argument("--annotated", type=Path, required=True, help="annotated_data.parquet (raw DR/DL for the curated panel)")
    parser.add_argument("--output-stem", type=Path, required=True, help="Output stem; pages land as <stem>_k<k>.pdf")
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name for the figure title")
    parser.add_argument("--reference-k", type=int, default=None, help="Add the hand-curated panel at this k")
    parser.add_argument("--curated", type=Path, default=None, help="manual_curated_clusters.tsv (Systematic ID, cluster)")
    parser.add_argument("--col-num", type=int, default=3, help="Subplot columns per page (default 3)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render pages, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotGridConfig(
            labels=args.labels,
            scaled=args.scaled,
            annotated=args.annotated,
            output_stem=args.output_stem,
            dataset=args.dataset,
            reference_k=args.reference_k,
            curated=args.curated,
            col_num=args.col_num,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
