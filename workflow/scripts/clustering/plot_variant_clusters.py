#!/usr/bin/env python3

"""
Plot Variant Clusters — DR/DL scatter for one finalize variant
==================================================================

Visualizes one finalize variant's cluster assignments in the (DR, DL) feature
space: one page for the final `cluster` (1..final_n_clusters, WT = lowest DR),
plus — for merge-type variants (auto_merge / manual_merge) that keep a
`raw_cluster` column — a second page for the pre-merge intermediate clusters.
The pages use the house cnsplots style via
workflow/src/figure_render/cluster_scatter.py, with one colour per cluster id
across every variant.

Input
-----
- final_clusters.tsv (any finalize variant: direct / auto_merge / grid / manual_merge)

Output
------
- cluster_scatter.pdf: 1 page (direct/grid, no raw_cluster) or 2 pages
  (auto_merge/manual_merge — intermediate page first, then final)
- cluster_scatter.review.png: page-1 preview at the same measured bounds

Usage
-----
    python plot_variant_clusters.py \\
        --final-clusters results/3b_clustering/{dataset}/{variant}/final_clusters.tsv \\
        --output results/3b_clustering/{dataset}/{variant}/cluster_scatter.pdf \\
        --variant-label kmeans_merge9_auto

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-21
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

# 3. Third-party Imports
import matplotlib as mpl

# 2. Data Processing Imports
import pandas as pd

mpl.use("Agg")  # headless: this script only writes a PDF, never displays
from loguru import logger
from matplotlib.figure import Figure

# 4. Local Imports (relative path resolution)
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from figure_render.cluster_scatter import build_cluster_scatter_figure
from figures import save_pdf_pages
from logging_setup import setup_logger


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PlotVariantClustersConfig:
    """Inputs, output, and the figure title label for one variant's cluster scatter."""

    final_clusters: Path
    output: Path
    variant_label: str

    def validate(self) -> None:
        """Raise ValueError if the input is missing, then ensure the output dir exists."""
        if not self.final_clusters.exists():
            raise ValueError(f"Required input not found: {self.final_clusters}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
# =============================================================================
# CORE LOGIC
# =============================================================================
def _cluster_figures(df: pd.DataFrame, variant_label: str) -> list[Figure]:
    """Build the intermediate (when raw_cluster exists) and final cluster scatter pages."""
    figures: list[Figure] = []
    if "raw_cluster" in df.columns:
        intermediate = df.dropna(subset=["raw_cluster"]).copy()
        intermediate["raw_cluster"] = intermediate["raw_cluster"].astype(int)
        figure = build_cluster_scatter_figure(
            intermediate,
            cluster_col="raw_cluster",
            one_based=False,
            title=(
                f"{variant_label} — intermediate "
                f"(raw_cluster, n={intermediate['raw_cluster'].nunique()})"
            ),
        )
        figures.append(figure)

    final = df.dropna(subset=["cluster"]).copy()
    final["cluster"] = final["cluster"].astype(int)
    figure = build_cluster_scatter_figure(
        final,
        cluster_col="cluster",
        one_based=True,
        show_box=True,
        legend=True,
        title=f"{variant_label} — final (cluster, n={final['cluster'].nunique()})",
    )
    figures.append(figure)
    return figures


@logger.catch(reraise=True)
def run(config: PlotVariantClustersConfig) -> None:
    """Render the intermediate (if present) and final cluster scatters into one PDF."""
    config.validate()
    df = pd.read_csv(config.final_clusters, sep="\t", index_col=0)

    figures = _cluster_figures(df, config.variant_label)
    save_pdf_pages(config.output.with_suffix(""), figures)

    logger.success(
        f"Wrote {len(figures)}-page cluster scatter for {config.variant_label} -> {config.output}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Plot DR/DL cluster scatter for one finalize variant"
    )
    parser.add_argument(
        "--final-clusters",
        type=Path,
        required=True,
        help="final_clusters.tsv (any finalize variant)",
    )
    parser.add_argument("--output", type=Path, required=True, help="Output cluster_scatter.pdf")
    parser.add_argument(
        "--variant-label", type=str, required=True, help="Variant name, used in figure titles"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging"
    )
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render the figure, report the outcome."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotVariantClustersConfig(
            final_clusters=args.final_clusters,
            output=args.output,
            variant_label=args.variant_label,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
