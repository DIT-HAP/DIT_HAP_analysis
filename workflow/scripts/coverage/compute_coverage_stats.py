#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gene Insertion Coverage Statistics
====================================

Per-dataset: computes insertion coverage (in-gene vs intergenic, using the
source notebook's exact IN_GENE_FILTER quirk) and gene coverage (covered vs
not covered by DR, split by essentiality) from the release/ curve-fitting
outputs. Ported from
DIT_HAP_pipeline/workflow/notebooks/gene_coverage_analysis.ipynb.

Input
-----
- Insertion-level fitting_results.tsv (MultiIndex [Chr, Coordinate, Strand,
  Target]) — defines the total insertion set.
- Insertion-level annotations.tsv(.gz) (same MultiIndex, plus Type /
  Distance_to_stop_codon / Systematic ID) — carries the in-gene/intergenic
  call per insertion. NOTE: this table can have duplicate index entries
  (multiple annotated Features per coordinate, e.g. CDS + overlapping
  intron); duplicates are collapsed (any Feature passing IN_GENE_FILTER
  wins) before joining against fitting_results, so counts are byte-faithful
  to the notebook's `fitting_results.index.isin(annotations.query(...).index)`
  approach without inflating per-chromosome or per-gene counts.
- Gene-level fitting_results.tsv (Systematic ID, DR, DeletionLibrary_essentiality
  already native columns — no extra essentiality merge needed here, unlike
  clustering.smk's RevisedDeletion_essentiality injection). Legacy releases
  still ship the pre-rename um/lam headers instead of DR/DL; normalized on
  load (same quirk as workflow/src/clustering/candidates.py).

Output
------
- coverage_stats.tsv: one row per coverage metric (insertion/gene/essential/
  non_essential x total/covered-or-in_gene/not_covered-or-intergenic).
- coverage_figures.pdf: donut charts (insertion + gene + essential +
  non-essential coverage) + per-chromosome insertion coverage bars + DR/DL
  histograms split by essentiality (3 rows x 2 cols).

Usage
-----
    python compute_coverage_stats.py \\
        --fitting-results .../release/insertion_level/fitting_results.tsv \\
        --annotations .../release/insertion_level/annotations.tsv.gz \\
        --gene-level .../release/gene_level/fitting_results.tsv \\
        --output-stats results/coverage/{dataset}/coverage_stats.tsv \\
        --output-figures results/coverage/{dataset}/coverage_figures.pdf

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-19
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
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.plotting.generic import donut_chart  # noqa: E402
from workflow.src.plotting.style import AX_HEIGHT, AX_WIDTH  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Byte-faithful to the source notebook's Config.in_gene_filter: an insertion
# counts as "in a gene" only if it's annotated as non-intergenic AND at least
# 5bp upstream of the stop codon (the >4 threshold, not >=5, is the notebook's
# own quirk — kept verbatim).
IN_GENE_FILTER = "Type != 'Intergenic region' and Distance_to_stop_codon > 4"

# Legacy -> current metric column names (same quirk as
# workflow/src/clustering/candidates.py's _LEGACY_METRIC_RENAME): some
# datasets' gene-level fitting_results.tsv still ship the pre-rename um/lam
# headers instead of DR/DL.
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}

# Donut chart colors, byte-faithful to the notebook's per-chart hardcoded values.
_INSERTION_COVERAGE_COLORS = ["#c4954b", "#C0C0C0"]
_GENE_COVERAGE_COLORS = ["#6b99df", "#C0C0C0"]
_ESSENTIAL_COVERAGE_COLORS = ["#dd8369", "#C0C0C0"]
_NON_ESSENTIAL_COVERAGE_COLORS = ["#98a64e", "#C0C0C0"]

# DR/DL histogram bin edges + x-limits + per-essentiality row colors, byte-faithful
# to the notebook's "DR DL Histogram" cell.
_DR_BINS = np.arange(-0.2, 1.5, 0.05)
_DR_XLIM = (-0.2, 1.5)
_DL_BINS = np.arange(0, 15, 0.5)
_DL_XLIM = (0, 15)
_HIST_ROW_COLORS = ["#6b99df", "#dd8369", "#98a64e"]
_HIST_ROW_QUERIES = [
    "DeletionLibrary_essentiality.notna()",
    "DeletionLibrary_essentiality == 'E'",
    "DeletionLibrary_essentiality == 'V'",
]
_HIST_ROW_LABELS = ["All genes", "Essential", "Non-essential"]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class CoverageConfig:
    """Inputs and outputs for the gene insertion coverage analysis."""
    fitting_results: Path
    annotations: Path
    gene_level: Path
    output_stats: Path
    output_figures: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.fitting_results, self.annotations, self.gene_level]:
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


def load_gene_level(gene_level_path: Path) -> pd.DataFrame:
    """Load gene-level fitting statistics, normalizing legacy um/lam -> DR/DL columns."""
    gene_result = pd.read_csv(gene_level_path, sep="\t")
    rename = {
        old: new
        for old, new in _LEGACY_METRIC_RENAME.items()
        if old in gene_result.columns and new not in gene_result.columns
    }
    if rename:
        logger.info(f"Normalizing legacy metric columns: {rename}")
        gene_result = gene_result.rename(columns=rename)
    return gene_result


def load_insertion_level(fitting_results_path: Path, annotations_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load insertion-level fitting results + annotations, both indexed by [Chr, Coordinate, Strand, Target].

    The annotations table can carry duplicate index entries (multiple
    Features per coordinate, e.g. an overlapping CDS + intron record).
    Duplicates are collapsed to one row per index value — preferring a row
    that passes IN_GENE_FILTER when any duplicate does — then reindexed onto
    fitting_results' index, matching the notebook's `.isin()` semantics
    without inflating counts from the raw many-to-one annotation rows.
    """
    fitting_results = pd.read_csv(fitting_results_path, sep="\t", index_col=[0, 1, 2, 3])
    annotations = pd.read_csv(annotations_path, sep="\t", index_col=[0, 1, 2, 3])

    if annotations.index.duplicated().any():
        n_dup = annotations.index.duplicated().sum()
        logger.info(f"Collapsing {n_dup} duplicate-indexed annotation rows (keep in-gene pass if any)")
        passes = annotations.eval(IN_GENE_FILTER)
        annotations = (
            annotations.assign(_passes=passes)
            .sort_values("_passes", ascending=False)
            .loc[lambda df: ~df.index.duplicated(keep="first")]
            .drop(columns="_passes")
        )

    annotations = annotations.reindex(fitting_results.index)
    return fitting_results, annotations


# =============================================================================
# CORE LOGIC — coverage computations (unit-tested)
# =============================================================================
def compute_insertion_coverage(annotation: pd.DataFrame) -> dict[str, int]:
    """Count in-gene vs intergenic insertions by the exact IN_GENE_FILTER quirk."""
    total = len(annotation)
    in_gene = len(annotation.query(IN_GENE_FILTER))
    return {"total": total, "in_gene": in_gene, "intergenic": total - in_gene}


def compute_gene_coverage(gene_result: pd.DataFrame) -> dict[str, int]:
    """Count genes covered (DR not NaN) vs not covered (DR is NaN)."""
    total = len(gene_result)
    covered = len(gene_result.query("DR.notna()"))
    return {"total": total, "covered": covered, "not_covered": total - covered}


def compute_essentiality_coverage(gene_result: pd.DataFrame) -> dict[str, dict[str, int]]:
    """Split compute_gene_coverage by DeletionLibrary_essentiality == 'E' vs not."""
    essential = gene_result[gene_result["DeletionLibrary_essentiality"] == "E"]
    non_essential = gene_result[gene_result["DeletionLibrary_essentiality"] != "E"]
    return {
        "essential": compute_gene_coverage(essential),
        "non_essential": compute_gene_coverage(non_essential),
    }


def compute_per_chromosome_insertion_coverage(annotation: pd.DataFrame) -> pd.DataFrame:
    """Per-chromosome in-gene/intergenic insertion counts (Chr is the 1st index level)."""
    rows = []
    for chrom, group in annotation.groupby(level="Chr"):
        counts = compute_insertion_coverage(group)
        rows.append({"Chr": chrom, **counts})
    return pd.DataFrame(rows).sort_values("Chr").reset_index(drop=True)


# =============================================================================
# STATS TABLE ASSEMBLY
# =============================================================================
def build_stats_table(
    insertion_coverage: dict[str, int],
    gene_coverage: dict[str, int],
    essentiality_coverage: dict[str, dict[str, int]],
    per_chromosome: pd.DataFrame,
) -> pd.DataFrame:
    """Flatten all coverage dicts into one long-form stats table."""
    rows = [
        {"metric": "insertion", "category": "all", "total": insertion_coverage["total"],
         "covered": insertion_coverage["in_gene"], "not_covered": insertion_coverage["intergenic"]},
        {"metric": "gene", "category": "all", "total": gene_coverage["total"],
         "covered": gene_coverage["covered"], "not_covered": gene_coverage["not_covered"]},
        {"metric": "gene", "category": "essential", "total": essentiality_coverage["essential"]["total"],
         "covered": essentiality_coverage["essential"]["covered"],
         "not_covered": essentiality_coverage["essential"]["not_covered"]},
        {"metric": "gene", "category": "non_essential", "total": essentiality_coverage["non_essential"]["total"],
         "covered": essentiality_coverage["non_essential"]["covered"],
         "not_covered": essentiality_coverage["non_essential"]["not_covered"]},
    ]
    for _, row in per_chromosome.iterrows():
        rows.append({
            "metric": "insertion", "category": f"chr_{row['Chr']}", "total": row["total"],
            "covered": row["in_gene"], "not_covered": row["intergenic"],
        })
    return pd.DataFrame(rows)


# =============================================================================
# PLOTTING
# =============================================================================
def plot_coverage_donuts(
    insertion_coverage: dict[str, int],
    gene_coverage: dict[str, int],
    essentiality_coverage: dict[str, dict[str, int]],
    per_chromosome: pd.DataFrame,
) -> plt.Figure:
    """Donut charts for insertion/gene/essential/non-essential coverage + per-chromosome bars."""
    fig, axes = plt.subplot_mosaic(
        [["A", "C", "E"], ["B", "D", "E"]],
        figsize=(AX_WIDTH * 3, AX_HEIGHT * 2),
    )

    donut_chart(
        values=[insertion_coverage["in_gene"], insertion_coverage["intergenic"]],
        labels=["In genes", "Intergenic regions"],
        colors=_INSERTION_COVERAGE_COLORS,
        center_text=f"Total\n{insertion_coverage['total']:,}\ninsertions",
        ax=axes["A"],
    )
    axes["A"].set_title("Insertions in coding genes")

    donut_chart(
        values=[gene_coverage["covered"], gene_coverage["not_covered"]],
        labels=["Covered", "Not covered"],
        colors=_GENE_COVERAGE_COLORS,
        center_text=f"Total\n{gene_coverage['total']:,}\ngenes",
        ax=axes["B"],
    )
    axes["B"].set_title("Gene coverage by insertions")

    essential = essentiality_coverage["essential"]
    donut_chart(
        values=[essential["covered"], essential["not_covered"]],
        labels=["Covered", "Not covered"],
        colors=_ESSENTIAL_COVERAGE_COLORS,
        center_text=f"Total\n{essential['total']:,}\nessential\ngenes",
        ax=axes["C"],
    )
    axes["C"].set_title("Essential gene\ncoverage by insertions")

    non_essential = essentiality_coverage["non_essential"]
    donut_chart(
        values=[non_essential["covered"], non_essential["not_covered"]],
        labels=["Covered", "Not covered"],
        colors=_NON_ESSENTIAL_COVERAGE_COLORS,
        center_text=f"Total\n{non_essential['total']:,}\nnon-essential\ngenes",
        ax=axes["D"],
    )
    axes["D"].set_title("Non-essential gene\ncoverage by insertions")

    ax = axes["E"]
    x = np.arange(len(per_chromosome))
    ax.bar(x, per_chromosome["in_gene"], label="In genes", color=_INSERTION_COVERAGE_COLORS[0])
    ax.bar(x, per_chromosome["intergenic"], bottom=per_chromosome["in_gene"],
           label="Intergenic regions", color=_INSERTION_COVERAGE_COLORS[1])
    ax.set_xticks(x)
    ax.set_xticklabels(per_chromosome["Chr"])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Number of insertions")
    ax.set_title("Per-chromosome insertion coverage")
    ax.legend()

    fig.tight_layout(h_pad=1, w_pad=1)
    return fig


def plot_dr_dl_histograms(gene_result: pd.DataFrame) -> plt.Figure:
    """3 rows (all/essential/non-essential) x 2 cols (DR/DL) histogram grid."""
    fig, axes = plt.subplots(3, 2, figsize=(AX_WIDTH * 2, AX_HEIGHT * 2))

    for col, col_feature in enumerate(["DR", "DL"]):
        bins, xlim = (_DR_BINS, _DR_XLIM) if col_feature == "DR" else (_DL_BINS, _DL_XLIM)
        for row, row_query in enumerate(_HIST_ROW_QUERIES):
            ax = axes[row, col]
            data = gene_result.query(row_query)[col_feature].dropna()
            ax.hist(data, bins=bins, rwidth=0.9, color=_HIST_ROW_COLORS[row])
            ax.set_xlim(xlim)
            if col == 0:
                ax.set_ylabel(f"{_HIST_ROW_LABELS[row]}\nNumber of genes")
            if row == 0:
                ax.set_title(f"{col_feature} distribution")

    fig.tight_layout()
    return fig


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: CoverageConfig) -> None:
    """Load -> compute coverage stats -> save TSV + figures."""
    config.validate()

    gene_result = load_gene_level(config.gene_level)
    fitting_results, annotations = load_insertion_level(config.fitting_results, config.annotations)

    insertion_coverage = compute_insertion_coverage(annotations)
    gene_coverage = compute_gene_coverage(gene_result)
    essentiality_coverage = compute_essentiality_coverage(gene_result)
    per_chromosome = compute_per_chromosome_insertion_coverage(annotations)

    stats_table = build_stats_table(insertion_coverage, gene_coverage, essentiality_coverage, per_chromosome)
    stats_table.to_csv(config.output_stats, sep="\t", index=False)

    fig_donuts = plot_coverage_donuts(insertion_coverage, gene_coverage, essentiality_coverage, per_chromosome)
    fig_hist = plot_dr_dl_histograms(gene_result)

    with PdfPages(config.output_figures) as pdf:
        pdf.savefig(fig_donuts, dpi=300, bbox_inches="tight")
        pdf.savefig(fig_hist, dpi=300, bbox_inches="tight")
    plt.close(fig_donuts)
    plt.close(fig_hist)

    logger.success(
        f"Coverage: {insertion_coverage['in_gene']:,}/{insertion_coverage['total']:,} insertions in-gene, "
        f"{gene_coverage['covered']:,}/{gene_coverage['total']:,} genes covered "
        f"({essentiality_coverage['essential']['covered']:,}/{essentiality_coverage['essential']['total']:,} essential)"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Compute gene insertion coverage statistics")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Insertion-level fitting_results.tsv")
    parser.add_argument("--annotations", type=Path, required=True, help="Insertion-level annotations.tsv(.gz)")
    parser.add_argument("--gene-level", type=Path, required=True, help="Gene-level fitting_results.tsv")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output coverage stats TSV")
    parser.add_argument("--output-figures", type=Path, required=True, help="Output coverage figures PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = CoverageConfig(
            fitting_results=args.fitting_results,
            annotations=args.annotations,
            gene_level=args.gene_level,
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
