#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot All Variant Clusters — one-page grid of every variant's final clusters
=============================================================================

Companion to plot_variant_clusters.py (which writes ONE pdf per variant): this
gathers EVERY variant's final_clusters.tsv into a single grid PDF so the final
`cluster` assignments of all finalize strategies can be eyeballed side by side
on identical DR/DL axes. Each subplot reuses the shared per-cluster color mapping
(workflow.src.plotting.gene_level.plot_cluster_on_axis), so cluster colors are
comparable across panels and match the per-variant scatters.

Input
-----
- final_clusters.tsv for each variant (paths + labels passed in parallel lists)

Output
------
- all_variants_cluster_scatter.pdf: one page, a grid of DR/DL scatters
  (one subplot per variant, final `cluster` column, WT = lowest DR)

Usage
-----
    python plot_all_variant_clusters.py \\
        --final-clusters path/to/v1/final_clusters.tsv path/to/v2/final_clusters.tsv \\
        --variant-labels kmeans_direct9 kmeans_merge9_auto \\
        --output results/clustering/final/{dataset}/all_variants_cluster_scatter.pdf \\
        --dataset HD_DIT_HAP

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-22
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
import matplotlib

matplotlib.use("Agg")  # headless: this script only writes a PDF, never displays
import matplotlib.pyplot as plt  # noqa: E402
from loguru import logger  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.plotting.gene_level import plot_cluster_on_axis  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PlotAllVariantClustersConfig:
    """Parallel lists of variant tsv paths + labels, the dataset name, and the output PDF."""
    final_clusters: list[Path]
    variant_labels: list[str]
    dataset: str
    output: Path
    col_num: int = 3

    def validate(self) -> None:
        """Raise ValueError if lists mismatch or any input is missing, then ensure output dir exists."""
        if len(self.final_clusters) != len(self.variant_labels):
            raise ValueError(
                f"Got {len(self.final_clusters)} tsv paths but {len(self.variant_labels)} labels; must match"
            )
        missing = [str(p) for p in self.final_clusters if not p.exists()]
        if missing:
            raise ValueError(f"Required input(s) not found: {missing}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotAllVariantClustersConfig) -> None:
    """Render every variant's final `cluster` scatter into one grid page."""
    config.validate()

    n = len(config.variant_labels)
    col_num = min(config.col_num, n)
    row_num = int(np.ceil(n / col_num))

    fig, axes = plt.subplots(row_num, col_num, figsize=(6 * col_num, 5.5 * row_num), sharex=True, sharey=True)
    axes = np.atleast_1d(axes).flatten()

    for idx, (path, label) in enumerate(zip(config.final_clusters, config.variant_labels)):
        df = pd.read_csv(path, sep="\t", index_col=0)
        final = df.dropna(subset=["cluster"]).copy()
        final["cluster"] = final["cluster"].astype(int)
        plot_cluster_on_axis(
            axes[idx], final, "cluster", cluster_minus_one=True, show_box=True,
            title=f"{label} (n_clusters={final['cluster'].nunique()})",
        )
        logger.info(f"  {label}: {len(final)} genes, {final['cluster'].nunique()} clusters")

    # Hide any unused axes in the last row.
    for j in range(n, len(axes)):
        fig.delaxes(axes[j])

    fig.suptitle(f"{config.dataset} — final clusters across all variants", fontsize=16, y=1.0)
    plt.tight_layout()

    with PdfPages(config.output) as pdf:
        pdf.savefig(fig)
    plt.close(fig)

    logger.success(f"Wrote {n}-variant cluster-scatter grid -> {config.output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Plot every variant's final clusters into one grid PDF")
    parser.add_argument("--final-clusters", type=Path, nargs="+", required=True, help="final_clusters.tsv for each variant")
    parser.add_argument("--variant-labels", type=str, nargs="+", required=True, help="Variant names (parallel to --final-clusters)")
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name (used in the figure title)")
    parser.add_argument("--output", type=Path, required=True, help="Output all_variants_cluster_scatter.pdf")
    parser.add_argument("--col-num", type=int, default=3, help="Number of subplot columns in the grid")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render the grid, report the outcome."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotAllVariantClustersConfig(
            final_clusters=args.final_clusters,
            variant_labels=args.variant_labels,
            dataset=args.dataset,
            output=args.output,
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
