#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Variant Clusters — DR/DL scatter for one finalize variant
==================================================================

Visualizes one finalize variant's cluster assignments in the (DR, DL) feature
space: one page for the final `cluster` (1..final_n_clusters, WT = lowest DR),
plus — for merge-type variants (auto_merge / manual_merge) that keep a
`raw_cluster` column — a second page for the pre-merge intermediate clusters.
Different clusters get different colors (reuses the notebook's cluster-review
plot: workflow/src/plotting/gene_level.visualize_cluster_on_feature_space), so
the pipeline plot and the manual-merge notebook's review plot look identical.

Input
-----
- final_clusters.tsv (any finalize variant: direct / auto_merge / grid / manual_merge)

Output
------
- cluster_scatter.pdf: 1 page (direct/grid, no raw_cluster) or 2 pages
  (auto_merge/manual_merge — intermediate page first, then final)

Usage
-----
    python plot_variant_clusters.py \\
        --final-clusters results/clustering/final/{dataset}/{variant}/final_clusters.tsv \\
        --output results/clustering/final/{dataset}/{variant}/cluster_scatter.pdf \\
        --variant-label kmeans_merge9_auto

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-21
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
import pandas as pd

# 3. Third-party Imports
import matplotlib

matplotlib.use("Agg")  # headless: this script only writes a PDF, never displays
from loguru import logger  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.plotting.gene_level import visualize_cluster_on_feature_space  # noqa: E402


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
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotVariantClustersConfig) -> None:
    """Render the intermediate (if present) and final cluster scatters into one PDF."""
    config.validate()
    df = pd.read_csv(config.final_clusters, sep="\t", index_col=0)

    pages = 0
    with PdfPages(config.output) as pdf:
        if "raw_cluster" in df.columns:
            intermediate = df.dropna(subset=["raw_cluster"]).copy()
            intermediate["raw_cluster"] = intermediate["raw_cluster"].astype(int)
            fig = visualize_cluster_on_feature_space(intermediate, "raw_cluster")
            fig.suptitle(f"{config.variant_label} — intermediate (raw_cluster, n={intermediate['raw_cluster'].nunique()})")
            pdf.savefig(fig)
            fig.clf()
            pages += 1

        final = df.dropna(subset=["cluster"]).copy()
        final["cluster"] = final["cluster"].astype(int)
        fig = visualize_cluster_on_feature_space(final, "cluster", cluster_minus_one=True, show_box=True, legend=True)
        fig.suptitle(f"{config.variant_label} — final (cluster, n={final['cluster'].nunique()})")
        pdf.savefig(fig)
        fig.clf()
        pages += 1

    logger.success(f"Wrote {pages}-page cluster scatter for {config.variant_label} -> {config.output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Plot DR/DL cluster scatter for one finalize variant")
    parser.add_argument("--final-clusters", type=Path, required=True, help="final_clusters.tsv (any finalize variant)")
    parser.add_argument("--output", type=Path, required=True, help="Output cluster_scatter.pdf")
    parser.add_argument("--variant-label", type=str, required=True, help="Variant name, used in figure titles")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
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
