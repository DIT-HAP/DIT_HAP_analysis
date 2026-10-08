#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Clusters Grid — one page per k, all methods side by side
=============================================================

Renders the (method x k) label table from compute_clustering_grid as a
multi-page PDF: page k is a grid of DR/DL scatters, one subplot per method,
all on shared axes so a reader can see how the five methods partition the
same space at the same cluster count (the visual companion to
k_selection.pdf's numeric QC). Cluster colors are one shared mapping across
every page (DR-ascending ids), so the same cluster reads as the same colour
everywhere.

Input
-----
- grid_labels.parquet (from compute_clustering_grid)
- scaled_data.parquet (the DR/DL coordinates to scatter)

Output
------
- clusters_grid.pdf + .review.png: one page per k in the table.

Usage
-----
    python plot_clusters_grid.py \\
        --labels results/clustering/{dataset}/_work/grid_labels.parquet \\
        --scaled results/clustering/{dataset}/_work/scaled_data.parquet \\
        --output results/clustering/{dataset}/clusters_grid.pdf \\
        --dataset HD_DIT_HAP

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-10-08
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
import matplotlib

matplotlib.use("Agg")  # headless: writes PDF/PNG only
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from loguru import logger  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from figures import apply_house_style, grid_axes, save_dual  # noqa: E402
from logging_setup import setup_logger  # noqa: E402
from plotting.gene_level import CLUSTER_COLORS  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotGridConfig:
    """Label table + coordinates in, one multi-page PDF out."""
    labels: Path
    scaled: Path
    output: Path
    dataset: str
    col_num: int = 3

    def validate(self) -> None:
        """Raise ValueError if any input is missing, then ensure the output dir exists."""
        for path in [self.labels, self.scaled]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
@logger.catch(reraise=True)
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


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotGridConfig) -> None:
    """Render every k's method grid into one multi-page PDF (+ review PNG of the last page)."""
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
    col_num = min(config.col_num, len(methods))
    row_num = -(-len(methods) // col_num)  # ceiling division

    with PdfPages(config.output) as pdf:
        for k in sorted(labels["k"].unique()):
            # One panel per method; suptitle needs headroom so it doesn't sit on
            # row 0's panel labels.
            axes = grid_axes(row_num, col_num, labels=None)
            fig = plt.gcf()
            fig.subplots_adjust(top=0.90, hspace=0.45, wspace=0.25)
            for method, ax in zip(methods, axes):
                df = labels[(labels["method"] == method) & (labels["k"] == k)]
                plot_df = scaled.loc[df["gene"]].assign(cluster=df["cluster"].to_numpy())
                scatter_clusters(ax, plot_df, colors)
                ax.set_title(f"{method} (k={k}, n={len(plot_df)})", fontsize=7)
            for ax in axes[len(methods):]:
                ax.set_visible(False)
            fig.suptitle(f"{config.dataset} — all methods at k={k}", fontsize=10)
            pdf.savefig(fig)
            save_dual(config.output.parent / f"{config.output.name}_k{k}")
            plt.close(fig)
        logger.success(f"Wrote {len(labels['k'].unique())} pages -> {config.output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Plot one page per k with all methods side by side")
    parser.add_argument("--labels", type=Path, required=True, help="grid_labels.parquet from compute_clustering_grid")
    parser.add_argument("--scaled", type=Path, required=True, help="scaled_data.parquet from the spine")
    parser.add_argument("--output", type=Path, required=True, help="Output clusters_grid.pdf")
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name for the figure title")
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
            output=args.output,
            dataset=args.dataset,
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
