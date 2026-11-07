#!/usr/bin/env python3

"""
Plot All Variant Clusters — one-page grid of every variant's final clusters
=============================================================================

Companion to plot_variant_clusters.py (which writes ONE pdf per variant): this
gathers EVERY variant's final_clusters.tsv into a single grid PDF so the final
`cluster` assignments of all finalize strategies can be eyeballed side by side
on identical DR/DL axes. Panels use the house cnsplots style via
workflow/src/figure_render/cluster_scatter.py, with one colour per cluster id
across every variant.

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
        --output results/3b_clustering/{dataset}/all_variants_cluster_scatter.pdf \\
        --dataset HD_DIT_HAP

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-22
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
import matplotlib.pyplot as plt

# 4. Local Imports (relative path resolution)
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from figure_render.cluster_scatter import build_cluster_grid_figure
from figures import save_dual
from logging_setup import setup_logger


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
        """Raise ValueError if lists mismatch or any input is missing, then ensure output dir
        exists."""
        if len(self.final_clusters) != len(self.variant_labels):
            raise ValueError(
                f"Got {len(self.final_clusters)} tsv paths but "
                f"{len(self.variant_labels)} labels; must match"
            )
        missing = [str(p) for p in self.final_clusters if not p.exists()]
        if missing:
            raise ValueError(f"Required input(s) not found: {missing}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotAllVariantClustersConfig) -> None:
    """Render every variant's final `cluster` scatter into one grid page."""
    config.validate()

    frames: dict[str, pd.DataFrame] = {}
    for path, label in zip(config.final_clusters, config.variant_labels, strict=False):
        df = pd.read_csv(path, sep="\t", index_col=0)
        final = df.dropna(subset=["cluster"]).copy()
        final["cluster"] = final["cluster"].astype(int)
        frames[label] = final
        logger.info(f"  {label}: {len(final)} genes, {final['cluster'].nunique()} clusters")

    figure = build_cluster_grid_figure(
        frames,
        col_num=config.col_num,
        title=f"{config.dataset} — final clusters across all variants",
    )
    save_dual(config.output.with_suffix(""))
    plt.close(figure)

    logger.success(f"Wrote {len(frames)}-variant cluster-scatter grid -> {config.output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Plot every variant's final clusters into one grid PDF"
    )
    parser.add_argument(
        "--final-clusters",
        type=Path,
        nargs="+",
        required=True,
        help="final_clusters.tsv for each variant",
    )
    parser.add_argument(
        "--variant-labels",
        type=str,
        nargs="+",
        required=True,
        help="Variant names (parallel to --final-clusters)",
    )
    parser.add_argument(
        "--dataset", type=str, required=True, help="Dataset name (used in the figure title)"
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="Output all_variants_cluster_scatter.pdf"
    )
    parser.add_argument(
        "--col-num", type=int, default=3, help="Number of subplot columns in the grid"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging"
    )
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
