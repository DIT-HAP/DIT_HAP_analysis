#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gene-Group Coherence Visualization
===================================

Reads computed coherence metrics from Parquet and generates a multi-panel
overview figure showing group-size distribution, z-score distribution, centroid
positions in fitness space, and coherence-vs-biology panels.

This is the plotting companion to compute_coherence.py, part of the
computation-plotting decoupling refactor (ADR-0001).

Input
-----
- coherence.parquet: computed metrics from compute_coherence.py, one row per
  group with columns: source, group_id, group_name, term_size, n_group_genes,
  covered_genes, centroid_x, centroid_y, median_distance, mean_distance,
  std_distance, min_distance, max_distance, mpd, z_score, p_value,
  n_permutations, mean_pairwise_distance_zscore, mean_pairwise_distance_p_value,
  p_fdr
- group_annotation_long.tsv: the prepared unified long-table (for panel A
  shared-subunit fraction)
- features.tsv: optional gene features table (for panels B/D biology metrics)

Output
------
- coherence.pdf: multi-panel figure with histograms, centroid map, and
  coherence-vs-biology scatter panels

Usage
-----
    python plot_coherence.py \\
        --input results/coherence/{dataset}/{source}/coherence.parquet \\
        --annotation results/coherence/{dataset}/{source}/group_annotation_long.tsv \\
        --output results/coherence/{dataset}/{source}/coherence.pdf \\
        [--features results/features/{pombase_version}/pombe_coding_gene_protein_features.tsv]

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-09-03
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Third-party
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from loguru import logger  # noqa: E402

# 3. Local
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402
from plotting.style import AX_HEIGHT, AX_WIDTH  # noqa: E402


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PlotConfig:
    """Inputs and outputs for coherence visualization."""
    input_metrics: Path
    annotation: Path
    output: Path
    features: Path | None = None

    def validate(self) -> None:
        """Raise ValueError if inputs missing, then create output dir."""
        required = [self.input_metrics, self.annotation]
        if self.features is not None:
            required.append(self.features)
        for path in required:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# BIOLOGY METRICS (from compute_coherence.py)
# =============================================================================
_ABUNDANCE_FEATURE_CANDIDATES = [
    "copies_per_cell_EMM_Proliferating_Cell",
    "copies_per_cell_EMMN_Quiescent_Cell",
    "mean_EMM_Proliferating_Cell_RNA_Abundance",
]
_CONSERVATION_FEATURE_CANDIDATES = ["evolutionary_rate"]


def _first_present_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first candidate column present in df.columns, or None."""
    for col in candidates:
        if col in df.columns:
            return col
    return None


def shared_subunit_fraction(long_table: pd.DataFrame) -> dict[str, float]:
    """Per group, fraction of members also in other groups (shared-subunit metric)."""
    member_to_groups = long_table.groupby("Systematic ID")["group_id"].apply(set).to_dict()
    result = {}
    for group_id, group_df in long_table.groupby("group_id"):
        members = group_df["Systematic ID"].unique()
        n_members = len(members)
        if n_members == 0:
            result[group_id] = 0.0
            continue
        shared = sum(1 for gene in members if len(member_to_groups.get(gene, set())) > 1)
        result[group_id] = shared / n_members
    return result


def member_feature_cv(long_table: pd.DataFrame, features: pd.DataFrame, feature_col: str) -> dict[str, float]:
    """Per group, coefficient of variation of a numeric feature across members."""
    features_dict = features.set_index("gene_systematic_id")[feature_col].to_dict()
    result = {}
    for group_id, group_df in long_table.groupby("group_id"):
        members = group_df["Systematic ID"].unique()
        values = [features_dict[gene] for gene in members if gene in features_dict]
        values = [v for v in values if pd.notna(v) and v > 0]
        if len(values) < 2:
            result[group_id] = np.nan
        else:
            result[group_id] = float(np.std(values, ddof=1) / np.mean(values))
    return result


# =============================================================================
# PLOTTING
# =============================================================================
def plot_coherence(
    table: pd.DataFrame,
    long_table: pd.DataFrame,
    features: pd.DataFrame | None = None
) -> plt.Figure:
    """Multi-panel coherence overview: size/z histograms, centroid map, coherence-vs-biology panels.

    Panels: (1) term-size histogram + z-score histogram; (2) centroid-position
    scatter (x=centroid_x/"typical DR", y=centroid_y/"typical DL", color=z-score,
    size proportional to term_size); (3) Panel A, coherence vs shared-subunit
    fraction (always drawn); and — only when `features` is provided — (4) Panel B,
    coherence vs protein/RNA-abundance uniformity, and Panel D, coherence vs
    conservation uniformity. Empty table keeps a graceful placeholder.
    """
    if table.empty:
        fig, ax = plt.subplots(1, 1, figsize=(AX_WIDTH, AX_HEIGHT))
        ax.text(0.5, 0.5, "No groups passed the size filter", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        return fig

    # Panel A: always computable from the long-table alone.
    shared_frac = shared_subunit_fraction(long_table)

    # Panels B (abundance) and D (conservation): only when features are supplied
    biology_panels: list[tuple[str, dict[str, float], str, str]] = [
        ("A", shared_frac, "Shared-subunit fraction", "Coherence vs shared subunits"),
    ]
    if features is None:
        logger.warning("No features table provided; skipping biology panels B (abundance) and D (conservation)")
    if features is not None:
        abundance_col = _first_present_column(features, _ABUNDANCE_FEATURE_CANDIDATES)
        if abundance_col is None:
            logger.warning(
                f"No abundance feature column found (tried {_ABUNDANCE_FEATURE_CANDIDATES}); skipping panel B"
            )
        else:
            biology_panels.append(
                ("B", member_feature_cv(long_table, features, abundance_col),
                 f"Abundance CV ({abundance_col})", "Coherence vs abundance uniformity")
            )
        conservation_col = _first_present_column(features, _CONSERVATION_FEATURE_CANDIDATES)
        if conservation_col is None:
            logger.warning(
                f"No conservation feature column found (tried {_CONSERVATION_FEATURE_CANDIDATES}); skipping panel D"
            )
        else:
            biology_panels.append(
                ("D", member_feature_cv(long_table, features, conservation_col),
                 f"Conservation CV ({conservation_col})", "Coherence vs conservation uniformity")
            )

    # Total axes: size hist + z hist + centroid map + one per biology panel.
    n_panels = 3 + len(biology_panels)
    ncols = 3
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(AX_WIDTH * ncols, AX_HEIGHT * nrows), squeeze=False
    )
    flat = axes.ravel()

    # Panel 1a: term-size histogram with log-spaced bins.
    ax_size = flat[0]
    log_min = np.log10(table["term_size"].min())
    log_max = np.log10(table["term_size"].max())
    log_bins = np.logspace(log_min, log_max, 21)
    ax_size.hist(table["term_size"], bins=log_bins, rwidth=0.9, color="#6b99df")
    ax_size.set_xlabel("Group size (DR>threshold members)")
    ax_size.set_ylabel("Number of groups")
    ax_size.set_title("Group size distribution")
    ax_size.set_xscale("log")

    # Panel 1b: z-score histogram.
    ax_z = flat[1]
    ax_z.hist(table["z_score"], bins=20, rwidth=0.9, color="#6b99df")
    ax_z.axvline(0.0, color="gray", linestyle="--", linewidth=0.8)
    ax_z.set_xlabel("MPD z-score (negative = coherent)")
    ax_z.set_ylabel("Number of groups")
    ax_z.set_title("Coherence z-score distribution")

    # Panel 2: centroid-position map in normalized fitness space.
    ax_centroid = flat[2]
    sizes = np.log1p(table["term_size"].to_numpy(dtype=float)) * 20
    scatter = ax_centroid.scatter(
        table["centroid_x"], table["centroid_y"], c=table["z_score"], s=sizes,
        cmap="coolwarm_r", alpha=0.8, edgecolors="none",
    )
    ax_centroid.set_xlabel("typical DR")
    ax_centroid.set_ylabel("typical DL/10")
    ax_centroid.set_title("Group centroid positions")
    fig.colorbar(scatter, ax=ax_centroid, label="z-score")

    # Size legend: show representative group sizes.
    legend_sizes = [3, 10, 30, 100]
    legend_handles = [
        plt.scatter([], [], s=np.log1p(float(s)) * 20, color="gray", alpha=0.6, edgecolors="none")
        for s in legend_sizes
    ]
    ax_centroid.legend(
        legend_handles, [str(s) for s in legend_sizes],
        title="Group size", loc="upper left", frameon=True, fontsize="small"
    )

    # Panel 3+: coherence-vs-biology scatter plots.
    for i, (label, cv_map, x_label, title) in enumerate(biology_panels, start=3):
        ax = flat[i]
        x_values = [cv_map.get(gid, np.nan) for gid in table["group_id"]]
        y_values = table["z_score"].to_numpy()
        ax.scatter(x_values, y_values, s=10, alpha=0.6, edgecolors="none", color="#6b99df")
        ax.axhline(0.0, color="gray", linestyle="--", linewidth=0.8)
        ax.set_xlabel(x_label)
        ax.set_ylabel("z-score")
        ax.set_title(title)
        ax.text(0.02, 0.98, label, transform=ax.transAxes, fontsize=14, fontweight="bold",
                va="top", ha="left")

    # Remove unused axes.
    for j in range(n_panels, len(flat)):
        fig.delaxes(flat[j])

    fig.tight_layout()
    return fig


# =============================================================================
# MAIN EXECUTION
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotConfig) -> None:
    """Generate coherence visualization from computed metrics."""
    config.validate()

    # Load inputs
    table = read_parquet(config.input_metrics)
    long_table = pd.read_csv(config.annotation, sep="\t")
    features = pd.read_csv(config.features, sep="\t") if config.features else None

    logger.info(f"Loaded {len(table):,} groups from {config.input_metrics}")

    # Generate figure
    fig = plot_coherence(table, long_table, features)

    # Save as PDF
    with PdfPages(config.output) as pdf:
        pdf.savefig(fig, dpi=300, bbox_inches="tight")
    plt.close(fig)

    logger.success(f"Wrote coherence figure to {config.output}")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Generate coherence visualization")
    parser.add_argument("--input", dest="input_metrics", type=Path, required=True,
                        help="Input coherence metrics Parquet from compute_coherence.py")
    parser.add_argument("--annotation", type=Path, required=True,
                        help="Prepared group_annotation_long.tsv (for panel A)")
    parser.add_argument("--features", type=Path, default=None,
                        help="Optional gene features TSV (for panels B/D)")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output coherence figure PDF")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run plotting, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotConfig(
            input_metrics=args.input_metrics,
            annotation=args.annotation,
            output=args.output,
            features=args.features,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
