"""
Gene Insertion Coverage — Core Logic
=====================================

Shared constants, loaders, coverage computations, stats-table assembly, and
figure builders for the coverage stage. Ported from
DIT_HAP_pipeline/workflow/notebooks/gene_coverage_analysis.ipynb and factored
out of the original single-script port so the stage can be split into
independent Snakemake rules (prepare -> compute stats / plot figures), each
re-runnable on its own.

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

Usage
-----
    from workflow.src.coverage.core import (
        load_gene_level, load_insertion_level, resolve_duplicate_annotations,
        compute_insertion_coverage, compute_gene_coverage,
        compute_essentiality_coverage, compute_per_chromosome_insertion_coverage,
        build_stats_table, plot_coverage_donuts, plot_dr_dl_histograms,
    )
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
import matplotlib

matplotlib.use("Agg")  # headless: builders only write PDFs, never display
import matplotlib.pyplot as plt  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
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
# Essential/non-essential rows use the SAME == 'E' / == 'V' definition as
# compute_essentiality_coverage (see that function's docstring) — genes with
# DeletionLibrary_essentiality == 'Not_determined' land in neither row, only
# in the "All genes" (.notna()) row. Keep these two definitions in sync: a
# mismatch here previously caused coverage_stats.tsv and coverage_figures.pdf
# to report different non_essential totals for the same run.
_HIST_ROW_QUERIES = [
    "DeletionLibrary_essentiality.notna()",
    "DeletionLibrary_essentiality == 'E'",
    "DeletionLibrary_essentiality == 'V'",
]
_HIST_ROW_LABELS = ["All genes", "Essential", "Non-essential"]


# =============================================================================
# LOADERS
# =============================================================================
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


def resolve_duplicate_annotations(annotations: pd.DataFrame) -> pd.DataFrame:
    """Collapse duplicate-indexed annotation rows to one row per index value.

    The insertion-level annotations table can carry duplicate index entries
    (multiple Features per coordinate, e.g. an overlapping CDS + intron
    record). Among duplicates sharing an index value, the row that passes
    IN_GENE_FILTER wins if any duplicate does (matching the notebook's
    `.isin()` semantics: an insertion counts as in-gene if ANY of its
    annotation rows qualifies). Uses an explicit `kind="stable"` sort so the
    tie-break is deterministic rather than relying on pandas' default
    quicksort (which does not guarantee a stable order for equal keys): when
    no duplicate passes (or the whole group has no duplicates), the first
    row in the original file order is kept.
    """
    if not annotations.index.duplicated().any():
        return annotations

    n_dup = annotations.index.duplicated().sum()
    logger.info(f"Collapsing {n_dup} duplicate-indexed annotation rows (keep in-gene pass if any)")
    passes = annotations.eval(IN_GENE_FILTER)
    return (
        annotations.assign(_passes=passes)
        .sort_values("_passes", ascending=False, kind="stable")
        .loc[lambda df: ~df.index.duplicated(keep="first")]
        .drop(columns="_passes")
    )


def load_insertion_level(fitting_results_path: Path, annotations_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load insertion-level fitting results + annotations, both indexed by [Chr, Coordinate, Strand, Target].

    Annotation duplicates (see resolve_duplicate_annotations) are collapsed
    before reindexing onto fitting_results' index, so counts are byte-faithful
    to the notebook's `fitting_results.index.isin(annotations.query(...).index)`
    approach without inflating counts from the raw many-to-one annotation rows.
    """
    fitting_results = pd.read_csv(fitting_results_path, sep="\t", index_col=[0, 1, 2, 3])
    annotations = pd.read_csv(annotations_path, sep="\t", index_col=[0, 1, 2, 3])

    annotations = resolve_duplicate_annotations(annotations)
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
    """Split compute_gene_coverage by DeletionLibrary_essentiality == 'E' vs == 'V'.

    Byte-faithful to the source notebook, which only ever tested
    `== 'E'` / `== 'V'` (never `!= 'E'`). Some releases carry a third value,
    `Not_determined` (e.g. 198/4513 genes in HD_DIT_HAP) — those genes are
    EXCLUDED from both buckets here (previously an earlier draft folded them
    into "non_essential" via `!= 'E'`, which silently diverged from the
    `_HIST_ROW_QUERIES` == 'V' filter used by plot_dr_dl_histograms and
    produced inconsistent totals between coverage_stats.tsv and the PDF).
    """
    essential = gene_result[gene_result["DeletionLibrary_essentiality"] == "E"]
    non_essential = gene_result[gene_result["DeletionLibrary_essentiality"] == "V"]
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
        # Some chromosome names already start with "chr_" (e.g.
        # "chr_II_telomeric_gap") — avoid doubling the prefix into
        # "chr_chr_II_telomeric_gap".
        chr_label = row["Chr"] if str(row["Chr"]).startswith("chr_") else f"chr_{row['Chr']}"
        rows.append({
            "metric": "insertion", "category": chr_label, "total": row["total"],
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
