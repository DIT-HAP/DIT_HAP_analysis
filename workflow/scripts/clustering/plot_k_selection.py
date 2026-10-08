#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot K Selection — cluster-count QC panels
===========================================

Renders the per-(method, k) QC metrics from compute_clustering_grid as one
panel per metric, all methods as separate lines against k on shared axes.
This figure is the "which k / which method" evidence: read it next to
clusters_grid.pdf, whose scatter pages show the partitions the numbers score.

Input
-----
- k_selection_metrics.tsv (from compute_clustering_grid)

Output
------
- k_selection.pdf + .review.png: one page, one panel per metric.

Usage
-----
    python plot_k_selection.py \\
        --metrics results/clustering/{dataset}/_work/k_selection_metrics.tsv \\
        --output results/clustering/{dataset}/k_selection.pdf \\
        --k-reference 9

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

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from figures import apply_house_style, grid_axes, save_dual  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# The metrics to panel, in reading order, with a human axis title each.
# (silhouette / CH / DB are the three the old variant metrics table carried;
# plotting them per method turns the single-build table into a comparison.)
PANELS: tuple[tuple[str, str], ...] = (
    ("silhouette_score", "Silhouette (↑ better)"),
    ("calinski_harabasz_score", "Calinski-Harabasz (↑ better)"),
    ("davies_bouldin_score", "Davies-Bouldin (↓ better)"),
    ("n_clusters", "Clusters found"),
)


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class KSelConfig:
    """Metrics table in, one QC page out."""
    metrics: Path
    output: Path
    k_reference: int | None = None
    col_num: int = 2

    def validate(self) -> None:
        """Raise ValueError if the input is missing or has unknown metric columns."""
        if not self.metrics.exists():
            raise ValueError(f"Required input not found: {self.metrics}")
        table = pd.read_csv(self.metrics, sep="\t", nrows=1)
        missing = [name for name, _ in PANELS if name not in table.columns]
        if missing:
            raise ValueError(f"Metrics table lacks columns: {missing}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: KSelConfig) -> None:
    """Render one panel per QC metric, one line per method, against k."""
    config.validate()
    metrics = pd.read_csv(config.metrics, sep="\t")
    methods = sorted(metrics["method"].unique())
    n_panels = len(PANELS)
    col_num = min(config.col_num, n_panels)
    row_num = -(-n_panels // col_num)

    apply_house_style()
    axes = grid_axes(row_num, col_num)
    fig = plt.gcf()
    fig.subplots_adjust(hspace=0.45, wspace=0.35)

    for (column, title), ax in zip(PANELS, axes):
        for method in methods:
            subset = metrics[metrics["method"] == method].sort_values("k")
            ax.plot(subset["k"], subset[column], marker="o", markersize=2.5, label=method)
        ax.set_title(title)
        ax.set_xlabel("k")
        ax.set_ylabel("")
        if config.k_reference is not None:
            ax.axvline(config.k_reference, color="grey", linestyle=":", linewidth=0.8)

    handles, names = axes[0].get_legend_handles_labels()
    fig.legend(handles, names, loc="outside upper right", ncols=len(methods), fontsize=6)
    save_dual(config.output.with_suffix(""))
    logger.success(
        f"k_selection -> {config.output} ({n_panels} metric panels, {len(methods)} methods, "
        f"reference k={config.k_reference})"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Plot cluster-count QC panels per method")
    parser.add_argument("--metrics", type=Path, required=True, help="k_selection_metrics.tsv from compute_clustering_grid")
    parser.add_argument("--output", type=Path, required=True, help="Output k_selection.pdf")
    parser.add_argument("--k-reference", type=int, default=None, help="Mark this k with a vertical guide line")
    parser.add_argument("--col-num", type=int, default=2, help="Panel columns (default 2)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render panels, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = KSelConfig(metrics=args.metrics, output=args.output, k_reference=args.k_reference, col_num=args.col_num)
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
