#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Cluster Sizes — how many genes per cluster, per method and k
=================================================================

Renders the per-(method, k, cluster) gene counts from compute_clustering_grid's
profile table as one figure: a matrix per method (rows = k, columns = cluster
id, colour/annotated value = gene count), so a reader can see at a glance which
methods split the gene set into comparable blocks and which spawn slivers
(1-2 gene clusters) that the internal metrics average away. Also writes the
counts as a tidy TSV beside the figure.

Input
-----
- k_selection_profile.tsv (from compute_clustering_grid)

Output
------
- cluster_sizes.tsv: one row per (method, k, cluster) with n_genes / share
- cluster_sizes.pdf + .review.png (via save_dual): one panel per method

Usage
-----
    python plot_cluster_sizes.py \\
        --profile results/clustering/{dataset}/_work/k_selection_profile.tsv \\
        --output-table results/clustering/{dataset}/_work/cluster_sizes.tsv \\
        --output-figure results/clustering/{dataset}/cluster_sizes.pdf \\
        --col-num 3

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

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from figures import apply_house_style, grid_axes, panel_labels, save_dual  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class SizesConfig:
    """Profile table in; a tidy counts TSV + one panel-per-method figure out."""
    profile: Path
    output_table: Path
    output_figure: Path
    col_num: int = 3

    def validate(self) -> None:
        """Raise ValueError if the input is missing or lacks the count column, then ensure output dirs exist."""
        if not self.profile.exists():
            raise ValueError(f"Required input not found: {self.profile}")
        header = pd.read_csv(self.profile, sep="\t", nrows=1)
        if "n_genes" not in header.columns:
            raise ValueError(f"Profile table lacks the n_genes column: {header.columns.tolist()}")
        self.output_table.parent.mkdir(parents=True, exist_ok=True)
        self.output_figure.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: SizesConfig) -> None:
    """Write the tidy counts table and render one count-matrix panel per method."""
    config.validate()

    profile = pd.read_csv(config.profile, sep="\t")
    sizes = profile[["method", "k", "cluster", "n_genes", "share"]]
    sizes.sort_values(["method", "k", "cluster"]).to_csv(config.output_table, sep="\t", index=False)
    methods = sorted(profile["method"].unique())
    max_k = int(profile["k"].max())

    apply_house_style()
    col_num = min(config.col_num, len(methods))
    row_num = -(-len(methods) // col_num)
    labels = panel_labels(len(methods))
    axes = grid_axes(row_num, col_num, labels=labels)
    fig = plt.gcf()
    fig.subplots_adjust(hspace=0.55, wspace=0.30)

    for method, ax in zip(methods, axes):
        matrix = (
            profile[profile["method"] == method]
            .pivot_table(index="k", columns="cluster", values="n_genes")
            .reindex(columns=range(1, max_k + 1))
        )
        # imshow with a linear scale: counts span ~1..3000, and the slivers the
        # figure exists to expose are the LOW end — linear keeps them visible
        # where a log ramp would flatten everything above ~10 together.
        image = ax.imshow(matrix.to_numpy(), aspect="auto", cmap="viridis", vmin=0)
        ax.set_xticks(range(max_k), [str(c) for c in range(1, max_k + 1)], fontsize=5)
        ax.set_yticks(range(len(matrix.index)), [str(k) for k in matrix.index], fontsize=5)
        ax.set_title(f"{method}", fontsize=7)
        ax.set_xlabel("cluster id (DR-ascending)")
        ax.set_ylabel("k")

    # One shared colorbar in its own inset on the last (unused) grid cell —
    # anchoring to the panel axes list resizes every panel past fixing.
    cax = axes[-1]
    cax.axis("off")
    cbar = fig.colorbar(image, ax=cax, fraction=0.9, shrink=0.85, pad=0.0)
    cbar.set_label("genes per cluster")

    save_dual(config.output_figure.with_suffix(""))
    logger.success(
        f"cluster sizes -> {config.output_table} + {config.output_figure} "
        f"({len(methods)} method panels, k up to {max_k})"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render per-method cluster-count matrices")
    parser.add_argument("--profile", type=Path, required=True, help="k_selection_profile.tsv from compute_clustering_grid")
    parser.add_argument("--output-table", type=Path, required=True, help="Output cluster_sizes.tsv")
    parser.add_argument("--output-figure", type=Path, required=True, help="Output cluster_sizes.pdf")
    parser.add_argument("--col-num", type=int, default=3, help="Panel columns (default 3)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, write table, render figure."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = SizesConfig(
            profile=args.profile,
            output_table=args.output_table,
            output_figure=args.output_figure,
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
