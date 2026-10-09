#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Cluster Genealogy — how clusters relate from k=4 to k=15 per method
========================================================================

Replaces the flat count matrices of the first cluster-sizes figure: one
alluvial ("sankey") panel per method. Each column is a k; the bars in it are
that k's clusters, stacked by size; a band joins a cluster to the cluster at
the next k that holds most of its genes. Reading across, a bar that stays one
band survived the finer k; a bar that fans into several is where the extra
clusters came from — the question the flat matrix could only half answer.

Links above the minimum overlap are table-computed (see
clustering.genealogy.cluster_genealogy); the render is the same stacked
alluvial the v2 sweep shipped, ported onto the grid-search label table.

Input
-----
- grid_labels.parquet (from compute_clustering_grid)

Output
------
- cluster_genealogy.pdf + .review.png (via save_dual): one panel per method

Usage
-----
    python plot_cluster_sizes.py \\
        --labels results/clustering/{dataset}/_work/grid_labels.parquet \\
        --output results/clustering/{dataset}/cluster_genealogy.pdf \\
        --col-num 3

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-10-09
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

matplotlib.use("Agg")
import pandas as pd  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from clustering.genealogy import cluster_genealogy, render_genealogy_panel  # noqa: E402
from figures import PanelShape, apply_house_style, fit_panels, grid_axes, panel_labels, save_dual  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class GenealogyConfig:
    """Label table in; one alluvial panel per method out."""
    labels: Path
    output: Path
    col_num: int = 3

    def validate(self) -> None:
        """Raise ValueError if the label table is missing or malformed, then ensure the output dir exists."""
        if not self.labels.exists():
            raise ValueError(f"Required input not found: {self.labels}")
        schema = pd.read_parquet(self.labels, columns=["method", "k", "cluster"]).pipe(lambda frame: frame.head(0))
        missing = [name for name in ("method", "k", "cluster") if name not in schema.columns]
        if missing:
            raise ValueError(f"Label table lacks columns: {missing}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: GenealogyConfig) -> None:
    """Compute the per-method genealogy and render one alluvial panel per method."""
    config.validate()

    labels = pd.read_parquet(config.labels)
    k_values = sorted(labels["k"].unique())
    methods = sorted(labels["method"].unique())

    apply_house_style()
    col_num = min(config.col_num, len(methods))
    row_num = -(-len(methods) // col_num)
    axes = grid_axes(row_num, col_num, labels=panel_labels(len(methods)), shape=PanelShape.WIDE)

    for method, ax in zip(methods, axes):
        genealogy = cluster_genealogy(labels[labels["method"] == method], k_values)
        render_genealogy_panel(ax, genealogy, k_values=k_values, method=method)

    for ax in axes[len(methods):]:  # hide the unused grid cells
        ax.set_visible(False)

    fit_panels()
    save_dual(config.output.with_suffix(""))
    logger.success(f"cluster genealogy -> {config.output} ({len(methods)} method panels, k={k_values[0]}..{k_values[-1]})")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render per-method cluster genealogy alluvials")
    parser.add_argument("--labels", type=Path, required=True, help="grid_labels.parquet from compute_clustering_grid")
    parser.add_argument("--output", type=Path, required=True, help="Output cluster_genealogy.pdf")
    parser.add_argument("--col-num", type=int, default=3, help="Panel columns (default 3)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render panels, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = GenealogyConfig(labels=args.labels, output=args.output, col_num=args.col_num)
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
