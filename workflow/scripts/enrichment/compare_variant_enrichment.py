#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compare Variant Enrichment — Cross-variant enrichment comparison plots
=======================================================================

Compares GO enrichment results across multiple clustering finalize variants to
visualize how different clustering strategies affect functional enrichment. For
each cluster (1..final_n_clusters), plots the top N most significant GO terms
across all variants, showing p_fdr values as a heatmap and enrichment statistics
as bar charts.

Input
-----
- go_enrichment_full_filtered.tsv for each variant (from enrichment.smk)
- Variant names from config

Output
------
- variant_enrichment_comparison.pdf: multi-page comparison
  - Page 1: Top terms heatmap (cluster × variant)
  - Page 2: Enrichment counts per cluster per variant
  - Page 3+: Per-cluster detailed comparison

Usage
-----
    python compare_variant_enrichment.py \\
        --dataset HD_DIT_HAP \\
        --pombase-version 2026-06-01 \\
        --variants kmeans_direct9 kmeans_merge9_auto hc_direct9 \\
        --output results/enrichment/comparison/HD_DIT_HAP/variant_enrichment_comparison.pdf

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
import numpy as np

# 3. Third-party Imports
import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.plotting import style  # noqa: E402, F401


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class CompareEnrichmentConfig:
    """Inputs and output for cross-variant enrichment comparison."""
    dataset: str
    pombase_version: str
    variants: list[str]
    output: Path
    top_n_terms: int = 10
    fdr_threshold: float = 0.05

    def validate(self) -> None:
        """Ensure output directory exists."""
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


def load_variant_enrichment(dataset: str, variant: str, pombase_version: str) -> pd.DataFrame:
    """Load one variant's GO enrichment results and add variant column."""
    path = Path(f"results/enrichment/raw/{dataset}/{variant}/{pombase_version}/go_enrichment_full_filtered.tsv")
    if not path.exists():
        raise FileNotFoundError(f"Enrichment file not found: {path}")
    df = pd.read_csv(path, sep="\t")
    df["variant"] = variant
    return df


# =============================================================================
# PLOTTING
# =============================================================================
def plot_enrichment_heatmap(combined: pd.DataFrame, config: CompareEnrichmentConfig) -> plt.Figure:
    """Plot heatmap of -log10(p_fdr) for top terms across variants and clusters."""
    # Get top N most significant terms across all variants
    top_terms = (combined
                 .sort_values("p_fdr")
                 .drop_duplicates("term")
                 .head(config.top_n_terms)["term"]
                 .tolist())

    # Pivot: rows=terms, cols=variant×cluster
    subset = combined[combined["term"].isin(top_terms)].copy()
    subset["-log10_fdr"] = -np.log10(subset["p_fdr"].clip(lower=1e-50))
    subset["variant_cluster"] = subset["variant"] + "_C" + subset["Cluster"].astype(str)

    pivot = subset.pivot_table(
        index="term", columns="variant_cluster", values="-log10_fdr", fill_value=0
    )

    fig, ax = plt.subplots(figsize=(14, 8))
    sns.heatmap(pivot, cmap="YlOrRd", cbar_kws={"label": "-log10(FDR)"}, ax=ax, linewidths=0.5)
    ax.set_title(f"Top {config.top_n_terms} GO Terms Across Variants (HD_DIT_HAP)", fontsize=14, weight="bold")
    ax.set_xlabel("Variant × Cluster", fontsize=12)
    ax.set_ylabel("GO Term", fontsize=12)
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(fontsize=9)
    plt.tight_layout()
    return fig


def plot_enrichment_counts(combined: pd.DataFrame, config: CompareEnrichmentConfig) -> plt.Figure:
    """Plot number of significant terms per cluster per variant."""
    counts = (combined[combined["p_fdr"] < config.fdr_threshold]
              .groupby(["variant", "Cluster"])
              .size()
              .reset_index(name="count"))

    fig, ax = plt.subplots(figsize=(12, 6))
    pivot_counts = counts.pivot(index="Cluster", columns="variant", values="count").fillna(0)
    pivot_counts.plot(kind="bar", ax=ax, width=0.8)
    ax.set_title(f"Significant GO Terms per Cluster (FDR < {config.fdr_threshold})", fontsize=14, weight="bold")
    ax.set_xlabel("Cluster", fontsize=12)
    ax.set_ylabel("Number of Significant Terms", fontsize=12)
    ax.legend(title="Variant", fontsize=9)
    plt.xticks(rotation=0)
    plt.tight_layout()
    return fig


def plot_per_cluster_comparison(combined: pd.DataFrame, cluster: int, config: CompareEnrichmentConfig) -> plt.Figure:
    """Plot top N terms for one cluster across all variants."""
    cluster_data = combined[combined["Cluster"] == cluster].copy()

    # Get top N terms for this cluster
    top_terms = (cluster_data
                 .sort_values("p_fdr")
                 .drop_duplicates("term")
                 .head(config.top_n_terms)["term"]
                 .tolist())

    subset = cluster_data[cluster_data["term"].isin(top_terms)].copy()
    subset["-log10_fdr"] = -np.log10(subset["p_fdr"].clip(lower=1e-50))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Left: p_fdr comparison
    pivot_fdr = subset.pivot_table(index="term", columns="variant", values="-log10_fdr", fill_value=0)
    pivot_fdr.plot(kind="barh", ax=ax1, width=0.8)
    ax1.set_title(f"Cluster {cluster}: Top GO Terms by Significance", fontsize=12, weight="bold")
    ax1.set_xlabel("-log10(FDR)", fontsize=11)
    ax1.set_ylabel("")
    ax1.legend(title="Variant", fontsize=8, loc="lower right")

    # Right: gene ratio comparison
    pivot_ratio = subset.pivot_table(index="term", columns="variant", values="gene_ratio", fill_value=0)
    pivot_ratio.plot(kind="barh", ax=ax2, width=0.8)
    ax2.set_title(f"Cluster {cluster}: Gene Ratio", fontsize=12, weight="bold")
    ax2.set_xlabel("Gene Ratio", fontsize=11)
    ax2.set_ylabel("")
    ax2.legend(title="Variant", fontsize=8, loc="lower right")

    plt.tight_layout()
    return fig


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: CompareEnrichmentConfig) -> None:
    """Load enrichment results, generate comparison plots."""
    config.validate()

    # Load all variant enrichments
    logger.info(f"Loading enrichment results for {len(config.variants)} variants...")
    dfs = []
    for variant in config.variants:
        try:
            df = load_variant_enrichment(config.dataset, variant, config.pombase_version)
            dfs.append(df)
            logger.info(f"  {variant}: {len(df)} enriched terms")
        except FileNotFoundError as e:
            logger.warning(f"  {variant}: SKIPPED ({e})")

    if not dfs:
        raise ValueError("No enrichment files found for any variant")

    combined = pd.concat(dfs, ignore_index=True)
    clusters = sorted(combined["Cluster"].unique())
    logger.info(f"Combined {len(combined)} rows across {len(clusters)} clusters")

    # Generate plots
    pages = 0
    with PdfPages(config.output) as pdf:
        # Page 1: Overall heatmap
        fig = plot_enrichment_heatmap(combined, config)
        pdf.savefig(fig)
        plt.close(fig)
        pages += 1

        # Page 2: Counts per cluster
        fig = plot_enrichment_counts(combined, config)
        pdf.savefig(fig)
        plt.close(fig)
        pages += 1

        # Page 3+: Per-cluster comparisons
        for cluster in clusters:
            fig = plot_per_cluster_comparison(combined, cluster, config)
            pdf.savefig(fig)
            plt.close(fig)
            pages += 1

    logger.success(f"Wrote {pages}-page comparison to {config.output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Compare GO enrichment across clustering variants")
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name (e.g., HD_DIT_HAP)")
    parser.add_argument("--pombase-version", type=str, required=True, help="PomBase version (e.g., 2026-06-01)")
    parser.add_argument("--variants", nargs="+", required=True, help="Variant names to compare")
    parser.add_argument("--output", type=Path, required=True, help="Output PDF path")
    parser.add_argument("--top-n-terms", type=int, default=10, help="Number of top terms to show per cluster")
    parser.add_argument("--fdr-threshold", type=float, default=0.05, help="FDR significance threshold")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: load data, generate plots, report outcome."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = CompareEnrichmentConfig(
            dataset=args.dataset,
            pombase_version=args.pombase_version,
            variants=args.variants,
            output=args.output,
            top_n_terms=args.top_n_terms,
            fdr_threshold=args.fdr_threshold,
        )
        run(config)
    except (ValueError, FileNotFoundError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
