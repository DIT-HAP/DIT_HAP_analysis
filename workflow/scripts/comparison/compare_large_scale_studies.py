#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pairwise Fitness Comparison with Other Large-Scale Studies
============================================================

Per-dataset: merges the curated DIT-HAP gene clusters with the gRNA (HD data)
fitted parameters and the PomBase-derived protein-features table, then
correlates the DIT-HAP fitness metric against every other large-scale
fitness/depletion readout available (gRNA um, Barseq, transposon integration
density, colony size, growth rate, ...). Produces a Pearson r + p-value table
and a pairwise scatter matrix with a Gaussian-KDE density overlay per panel.
Ported from
DIT_HAP_pipeline/workflow/notebooks/compare_with_other_large_scale_studies.ipynb.

This is a simplified port: the notebook's altair repeat-grid, the KEGG BRITE
pathway jitter charts, and the per-GO-term feature-space PDFs are out of scope
here — only the clip-and-correlate fitness comparison + KDE scatter the task
calls for. The pairwise column selection is defensive: only fitness columns
actually present (with enough non-NaN data) are correlated/plotted, so the
script never KeyErrors on a schema that ships a subset of the study columns.

Input
-----
- final_clusters.tsv (Systematic ID + the DIT-HAP fitness metric; DR is the
  current name, legacy releases ship it as `um`, normalized on load) from the
  clustering finalize-variant stage, sourced via final_clusters_path(dataset,
  selected_variant). Only Systematic ID + DR are read here.
- gRNA HD-data fitted parameters TSV (Systematic ID + `um` gRNA fitness metric).
- pombe_coding_gene_protein_features.tsv (gene_systematic_id + the other
  large-scale study columns: Barseq_from_dulab/koch, integration density, ipkm,
  uipkm, colony_size_Malecki2016, Max Growth Rate, Colony Formation).

Output
------
- fitness_correlation_stats.tsv: long-form (col_x, col_y, pair, r, p_value, n),
  one row per unordered pair of available fitness columns.
- pairwise_fitness_comparison.pdf: scatter matrix (one panel per pair) with a
  Gaussian-KDE density overlay and Pearson r/p/n annotation.

Usage
-----
    python compare_large_scale_studies.py \\
        --final-clusters results/clustering/final/{dataset}/{variant}/final_clusters.tsv \\
        --protein-features results/features/{ver}/pombe_coding_gene_protein_features.tsv \\
        --grna-data resources/curated/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv \\
        --clip-upper 200 \\
        --output-stats results/comparison/{dataset}/fitness_correlation_stats.tsv \\
        --output-figures results/comparison/{dataset}/pairwise_fitness_comparison.pdf

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
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
from itertools import combinations
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, pearsonr

# 3. Third-party Imports
import matplotlib

matplotlib.use("Agg")  # headless: this script only writes a PDF, never displays
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.plotting.style import AX_HEIGHT, AX_WIDTH, COLORS  # noqa: E402

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Byte-faithful to the source notebook: the transposon integration-density
# metrics are heavy-tailed, so the notebook caps them at 200 with .clip(upper=200)
# before any plotting/correlation. Only these three columns are clipped.
CLIP_UPPER = 200
DENSITY_COLUMNS = [
    "Integration density, in-vivo (integrations/kb/million inserts)",
    "ipkm",
    "uipkm",
]

# Legacy -> current metric column names, same quirk as coverage.smk /
# verification / noncoding_rna: some curated final_clusters.tsv releases still
# ship the pre-rename `um`/`lam` headers instead of DR/DL.
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}

# The other large-scale study fitness/depletion columns to correlate against,
# byte-faithful to the notebook's fitness_data column list. These live on the
# protein-features table; the DIT-HAP metric (DR) and gRNA metric (um_gRNA) are
# merged in separately. Column selection at runtime is DEFENSIVE — only those
# actually present with enough non-NaN data are used (see select_fitness_columns).
STUDY_FITNESS_COLUMNS = [
    "Barseq_from_dulab",
    "Barseq_from_koch",
    "Integration density, in-vivo (integrations/kb/million inserts)",
    "ipkm",
    "uipkm",
    "colony_size_Malecki2016",
    "Max Growth Rate",
    "Colony Formation",
]

# The DIT-HAP and gRNA fitness metrics after merge (see build_fitness_table).
DIT_HAP_FITNESS_COLUMN = "um_DIT_HAP"
GRNA_FITNESS_COLUMN = "um_gRNA"

# A Pearson correlation needs at least this many complete (non-NaN) pairs to be
# meaningful; pairs below this are skipped (logged) rather than emitting a
# degenerate r/p that scipy warns or NaNs on.
MIN_PAIRS_FOR_CORRELATION = 3


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class ComparisonConfig:
    """Inputs, parameters, and outputs for the pairwise fitness comparison."""
    final_clusters: Path
    protein_features: Path
    grna_data: Path
    clip_upper: float
    output_stats: Path
    output_figures: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.final_clusters, self.protein_features, self.grna_data]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_stats, self.output_figures]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


def load_final_clusters(final_clusters_path: Path) -> pd.DataFrame:
    """Load the curated cluster table, normalizing legacy um/lam -> DR/DL columns."""
    clusters = pd.read_csv(final_clusters_path, sep="\t")
    rename = {
        old: new
        for old, new in _LEGACY_METRIC_RENAME.items()
        if old in clusters.columns and new not in clusters.columns
    }
    if rename:
        logger.info(f"Normalizing legacy metric columns: {rename}")
        clusters = clusters.rename(columns=rename)
    return clusters


# =============================================================================
# CORE LOGIC — clip + correlate (unit-tested)
# =============================================================================
def clip_density_columns(df: pd.DataFrame, clip_upper: float = CLIP_UPPER) -> pd.DataFrame:
    """Cap the three heavy-tailed integration-density columns at clip_upper.

    Byte-faithful to the notebook's .clip(upper=200) on exactly
    DENSITY_COLUMNS; every other column is left untouched. Columns absent from
    ``df`` are silently skipped so this works against a partial real schema.
    """
    result = df.copy()
    for column in DENSITY_COLUMNS:
        if column in result.columns:
            result[column] = result[column].clip(upper=clip_upper)
    return result


def compute_pearson_r(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    """Pearson (r, p_value) over the pairs where both x and y are non-NaN.

    Drops any pair with a NaN in either series before calling scipy's pearsonr
    (matches the notebook, which correlates only complete observations).
    """
    paired = pd.DataFrame({"x": x.to_numpy(), "y": y.to_numpy()}).dropna()
    r, p_value = pearsonr(paired["x"], paired["y"])
    return float(r), float(p_value)


# =============================================================================
# CORE LOGIC — merge + assembly
# =============================================================================
def build_fitness_table(
    final_clusters: pd.DataFrame,
    protein_features: pd.DataFrame,
    grna_data: pd.DataFrame,
    clip_upper: float = CLIP_UPPER,
) -> pd.DataFrame:
    """Merge protein-features + DIT-HAP clusters + gRNA into one fitness table.

    Mirrors the notebook's merged_fitness_data: the protein-features table
    (keyed on gene_systematic_id) is the spine, left-joined to the curated
    DIT-HAP metric and the gRNA metric on Systematic ID. Both metric columns
    are the fitting `DR`/`um` column, disambiguated to um_DIT_HAP / um_gRNA.
    The integration-density columns are clipped at clip_upper on the way out.
    """
    dit_hap_metric = _dit_hap_metric_column(final_clusters)
    grna_metric = _grna_metric_column(grna_data)

    merged = protein_features.merge(
        final_clusters[["Systematic ID", dit_hap_metric]].rename(
            columns={dit_hap_metric: DIT_HAP_FITNESS_COLUMN}
        ),
        left_on="gene_systematic_id",
        right_on="Systematic ID",
        how="left",
    ).merge(
        grna_data[["Systematic ID", grna_metric]].rename(
            columns={grna_metric: GRNA_FITNESS_COLUMN}
        ),
        left_on="gene_systematic_id",
        right_on="Systematic ID",
        how="left",
        suffixes=("_dithap", "_grna"),
    )
    return clip_density_columns(merged, clip_upper=clip_upper)


def _dit_hap_metric_column(final_clusters: pd.DataFrame) -> str:
    """Pick the DIT-HAP fitness metric column from the curated cluster table.

    Prefers the current `DR` name (load_final_clusters normalizes legacy `um`);
    falls back to `um` if a caller passed a raw frame. Raises if neither exists
    so a schema drift surfaces loudly instead of silently dropping the metric.
    """
    for candidate in ("DR", "um"):
        if candidate in final_clusters.columns:
            return candidate
    raise KeyError("final_clusters must contain a 'DR' (or legacy 'um') fitness column")


def _grna_metric_column(grna_data: pd.DataFrame) -> str:
    """Pick the gRNA fitness metric column (native `um`, or normalized `DR`)."""
    for candidate in ("um", "DR"):
        if candidate in grna_data.columns:
            return candidate
    raise KeyError("grna_data must contain a 'um' (or 'DR') fitness column")


def select_fitness_columns(fitness_table: pd.DataFrame) -> list[str]:
    """Return the fitness columns present with enough non-NaN data to correlate.

    Defensive against a partial real schema: only STUDY_FITNESS_COLUMNS plus the
    two merged metrics that actually exist in ``fitness_table`` and have at least
    MIN_PAIRS_FOR_CORRELATION non-NaN values are kept. Missing/too-sparse columns
    are logged and skipped so the pairwise loop can't KeyError at runtime.
    """
    candidates = STUDY_FITNESS_COLUMNS + [DIT_HAP_FITNESS_COLUMN, GRNA_FITNESS_COLUMN]
    available, missing = [], []
    for column in candidates:
        if column in fitness_table.columns and fitness_table[column].notna().sum() >= MIN_PAIRS_FOR_CORRELATION:
            available.append(column)
        else:
            missing.append(column)
    if missing:
        logger.warning(f"Skipping {len(missing)} fitness column(s) (absent or too sparse): {missing}")
    logger.info(f"Correlating {len(available)} fitness columns: {available}")
    return available


def compute_correlation_stats(fitness_table: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Long-form Pearson r/p_value/n for every unordered pair of ``columns``.

    Pairs with fewer than MIN_PAIRS_FOR_CORRELATION complete observations are
    skipped (logged) rather than emitting a degenerate r/p.
    """
    rows = []
    for col_x, col_y in combinations(columns, 2):
        paired = fitness_table[[col_x, col_y]].dropna()
        n = len(paired)
        if n < MIN_PAIRS_FOR_CORRELATION:
            logger.warning(f"Skipping pair ({col_x} vs {col_y}): only {n} complete pairs")
            continue
        r, p_value = compute_pearson_r(paired[col_x], paired[col_y])
        rows.append({
            "col_x": col_x,
            "col_y": col_y,
            "pair": f"{col_x} vs {col_y}",
            "r": r,
            "p_value": p_value,
            "n": n,
        })
    return pd.DataFrame(rows, columns=["col_x", "col_y", "pair", "r", "p_value", "n"])


# =============================================================================
# PLOTTING
# =============================================================================
def _plot_pair(ax: plt.Axes, fitness_table: pd.DataFrame, col_x: str, col_y: str) -> None:
    """Scatter of col_x vs col_y with a Gaussian-KDE density overlay + r/p/n text.

    Guarded by MIN_PAIRS_FOR_CORRELATION on the PAIRWISE overlap (not the
    per-column count select_fitness_columns uses): two columns can each be dense
    on their own yet share zero rows (e.g. real Barseq vs Max Growth Rate),
    which would make compute_pearson_r's scipy.pearsonr raise. Below-threshold
    panels are blanked instead — this also suppresses the misleading n=2 -> r=1
    panel. plot_pairwise_comparison already iterates only the surviving stats
    pairs, so this guard is belt-and-suspenders for direct callers.
    """
    paired = fitness_table[[col_x, col_y]].dropna()
    if len(paired) < MIN_PAIRS_FOR_CORRELATION:
        ax.axis("off")
        return
    x = paired[col_x].to_numpy()
    y = paired[col_y].to_numpy()

    ax.scatter(x, y, s=8, alpha=0.3, color=COLORS[0], edgecolors="none", zorder=1)

    # KDE contour overlay: needs >2 points and non-degenerate spread; a constant
    # column makes the covariance singular, so fall back to the bare scatter.
    if len(paired) > 2 and np.ptp(x) > 0 and np.ptp(y) > 0:
        try:
            xy = np.vstack([x, y])
            kde = gaussian_kde(xy)
            xi, yi = np.mgrid[x.min():x.max():60j, y.min():y.max():60j]
            zi = kde(np.vstack([xi.ravel(), yi.ravel()])).reshape(xi.shape)
            ax.contour(xi, yi, zi, levels=6, colors=COLORS[1], linewidths=0.6, zorder=2)
        except np.linalg.LinAlgError:
            logger.debug(f"KDE overlay skipped for ({col_x} vs {col_y}): singular covariance")

    r, p_value = compute_pearson_r(paired[col_x], paired[col_y])
    ax.set_title(f"r={r:.2f}, p={p_value:.1e}\nn={len(paired)}", fontsize=7)
    ax.set_xlabel(col_x, fontsize=6)
    ax.set_ylabel(col_y, fontsize=6)
    ax.tick_params(labelsize=6)


def plot_pairwise_comparison(fitness_table: pd.DataFrame, pairs: list[tuple[str, str]]) -> plt.Figure:
    """Scatter matrix (one panel per pair) with KDE overlay + r/p/n.

    ``pairs`` is the list of (col_x, col_y) that SURVIVED the per-pair overlap
    filter in compute_correlation_stats, so the PDF panel set always matches the
    stats TSV row set — the two outputs can't silently disagree. Grid is packed
    row-major so a non-square pair count leaves trailing axes blank (turned off)
    rather than skewing the layout.
    """
    n_pairs = max(len(pairs), 1)
    n_cols = min(4, n_pairs)
    n_rows = int(np.ceil(n_pairs / n_cols))

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(AX_WIDTH * n_cols, AX_HEIGHT * n_rows),
        squeeze=False,
    )
    flat_axes = axes.ravel()
    for ax, (col_x, col_y) in zip(flat_axes, pairs):
        _plot_pair(ax, fitness_table, col_x, col_y)
    for ax in flat_axes[len(pairs):]:
        ax.axis("off")

    fig.suptitle("Pairwise fitness comparison across large-scale studies", y=1.01)
    fig.tight_layout()
    return fig


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: ComparisonConfig) -> None:
    """Load -> merge -> clip -> correlate -> save TSV + scatter-matrix PDF."""
    config.validate()

    final_clusters = load_final_clusters(config.final_clusters)
    protein_features = pd.read_csv(config.protein_features, sep="\t")
    grna_data = pd.read_csv(config.grna_data, sep="\t")

    fitness_table = build_fitness_table(
        final_clusters, protein_features, grna_data, clip_upper=config.clip_upper
    )
    columns = select_fitness_columns(fitness_table)

    stats = compute_correlation_stats(fitness_table, columns)
    stats.to_csv(config.output_stats, sep="\t", index=False)

    # Drive the plot grid from the surviving stats pairs (post per-pair overlap
    # filter) so PDF panels and TSV rows always agree on which pairs exist.
    surviving_pairs = list(zip(stats["col_x"], stats["col_y"]))
    fig = plot_pairwise_comparison(fitness_table, surviving_pairs)
    with PdfPages(config.output_figures) as pdf:
        pdf.savefig(fig, dpi=300, bbox_inches="tight")
    plt.close(fig)

    logger.success(
        f"Comparison: {len(stats):,} fitness-column pairs correlated across "
        f"{len(columns):,} columns ({len(fitness_table):,} genes)"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Pairwise fitness comparison with other large-scale studies")
    parser.add_argument("--final-clusters", type=Path, required=True, help="Curated final_clusters.tsv")
    parser.add_argument("--protein-features", type=Path, required=True, help="pombe_coding_gene_protein_features.tsv")
    parser.add_argument("--grna-data", type=Path, required=True, help="gRNA HD-data fitted parameters TSV")
    parser.add_argument("--clip-upper", type=float, default=CLIP_UPPER, help="Upper cap for integration-density columns")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output correlation stats TSV")
    parser.add_argument("--output-figures", type=Path, required=True, help="Output pairwise comparison PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = ComparisonConfig(
            final_clusters=args.final_clusters,
            protein_features=args.protein_features,
            grna_data=args.grna_data,
            clip_upper=args.clip_upper,
            output_stats=args.output_stats,
            output_figures=args.output_figures,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())

