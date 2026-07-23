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
- Gene-level fitting_results.tsv (Systematic ID, DR, ... ). Legacy releases
  still ship the pre-rename um/lam headers instead of DR/DL; normalized on
  load (same quirk as workflow/src/clustering/candidates.py). Its native
  FYPOviability/DeletionLibrary_essentiality columns are dropped by
  prepare_coverage_data.py in favor of deletion_viability (gene_metadata) and
  essentiality (deletion_library_categories.xlsx) — same underlying facts,
  sourced for the FULL protein-coding gene universe instead of just the
  DIT-HAP-covered subset (see prepare_coverage_data.py's run() for the
  byte-for-byte equivalence check).

Usage
-----
    from workflow.src.coverage.core import (
        load_gene_level, load_insertion_level, resolve_duplicate_annotations,
        compute_insertion_coverage, compute_gene_coverage,
        compute_essentiality_coverage, compute_per_chromosome_insertion_coverage,
        compute_characterisation_status_coverage,
        compute_deletion_viability_coverage, compute_essentiality_category_coverage,
        build_stats_table, coverage_dicts_from_stats_table,
        plot_coverage_donuts, plot_dr_dl_histograms,
        plot_characterisation_status_donuts, plot_characterisation_status_histograms,
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
# characterisation_status donuts reuse the shared covered/not-covered scheme
# (covered = the gene-coverage blue, not-covered = the same grey as every other donut).
_CHARACTERISATION_COVERAGE_COLORS = ["#6b99df", "#C0C0C0"]

# DR/DL histogram bin edges + x-limits + per-essentiality row colors, byte-faithful
# to the notebook's "DR DL Histogram" cell.
_DR_BINS = np.arange(-0.2, 1.5, 0.05)
_DR_XLIM = (-0.2, 1.5)
_DL_BINS = np.arange(0, 15, 0.5)
_DL_XLIM = (0, 15)
_HIST_ROW_COLORS = ["#6b99df", "#dd8369", "#98a64e"]
# Essential/non-essential rows use the SAME == 'E' / == 'V' definition as
# compute_essentiality_coverage (see that function's docstring) — genes with
# essentiality == 'Not_determined' land in neither row, only in the
# "All genes" (.notna()) row. Keep these two definitions in sync: a mismatch
# here previously caused coverage_stats.tsv and coverage_figures.pdf to report
# different non_essential totals for the same run.
_HIST_ROW_QUERIES = [
    "essentiality.notna()",
    "essentiality == 'E'",
    "essentiality == 'V'",
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
# STATS-TABLE READBACK (so figures read the SAME numbers the stats rule wrote)
# =============================================================================
# Prefixes build_stats_table uses for per-category rows; readers strip them to
# recover the raw category label. Keep in sync with build_stats_table.
_CHARACTERISATION_PREFIX = "characterisation_"
_DELETION_VIABILITY_PREFIX = "deletion_viability_"
_ESSENTIALITY_CATEGORY_PREFIX = "essentiality_"


def coverage_dicts_from_stats_table(
    stats: pd.DataFrame,
) -> tuple[
    dict[str, int],
    dict[str, int],
    dict[str, dict[str, int]],
    pd.DataFrame,
    dict[str, dict[str, int]],
    dict[str, dict[str, int]],
    dict[str, dict[str, int]],
]:
    """Reconstruct the coverage dicts + per-chromosome table from a coverage_stats.tsv frame.

    Inverse of build_stats_table: lets plot_coverage_figures render donuts from the
    exact numbers compute_coverage_stats wrote, instead of recomputing them from the
    gene_result parquet (which risks figure/table drift if the two paths ever diverge).
    Returns (insertion_coverage, gene_coverage, essentiality_coverage, per_chromosome,
    characterisation_status_coverage, deletion_viability_coverage,
    essentiality_category_coverage) — the first four match plot_coverage_donuts' args
    (insertion dict re-exposes covered/not_covered as in_gene/intergenic), the
    characterisation_status dict feeds plot_characterisation_status_donuts, and the last
    two are the full per-category breakdowns (4 deletion_viability values, 3 essentiality
    values including Not_determined) that essentiality_coverage's essential/non_essential
    split omits.
    """
    def _row(metric: str, category: str) -> pd.Series:
        hit = stats[(stats["metric"] == metric) & (stats["category"] == category)]
        if hit.empty:
            raise ValueError(f"coverage_stats table missing required row: metric={metric!r}, category={category!r}")
        return hit.iloc[0]

    ins = _row("insertion", "all")
    insertion_coverage = {"total": int(ins["total"]), "in_gene": int(ins["covered"]), "intergenic": int(ins["not_covered"])}

    gene = _row("gene", "all")
    gene_coverage = {"total": int(gene["total"]), "covered": int(gene["covered"]), "not_covered": int(gene["not_covered"])}

    def _gene_cov(category: str) -> dict[str, int]:
        r = _row("gene", category)
        return {"total": int(r["total"]), "covered": int(r["covered"]), "not_covered": int(r["not_covered"])}

    essentiality_coverage = {"essential": _gene_cov("essential"), "non_essential": _gene_cov("non_essential")}

    # Per-chromosome insertion rows: any insertion row that isn't the "all" summary.
    per_chr_rows = stats[(stats["metric"] == "insertion") & (stats["category"] != "all")]
    per_chromosome = pd.DataFrame({
        "Chr": per_chr_rows["category"].str.replace(r"^chr_", "", regex=True),
        "total": per_chr_rows["total"].astype(int),
        "in_gene": per_chr_rows["covered"].astype(int),
        "intergenic": per_chr_rows["not_covered"].astype(int),
    }).reset_index(drop=True)

    def _prefixed_category_coverage(prefix: str) -> dict[str, dict[str, int]]:
        result: dict[str, dict[str, int]] = {}
        rows = stats[(stats["metric"] == "gene") & (stats["category"].str.startswith(prefix))]
        for _, r in rows.iterrows():
            key = r["category"][len(prefix):]
            result[key] = {"total": int(r["total"]), "covered": int(r["covered"]), "not_covered": int(r["not_covered"])}
        return result

    characterisation_status_coverage = _prefixed_category_coverage(_CHARACTERISATION_PREFIX)
    deletion_viability_coverage = _prefixed_category_coverage(_DELETION_VIABILITY_PREFIX)
    essentiality_category_coverage = _prefixed_category_coverage(_ESSENTIALITY_CATEGORY_PREFIX)

    return (
        insertion_coverage,
        gene_coverage,
        essentiality_coverage,
        per_chromosome,
        characterisation_status_coverage,
        deletion_viability_coverage,
        essentiality_category_coverage,
    )


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
    """Split compute_gene_coverage by essentiality == 'E' vs == 'V'.

    Byte-faithful to the source notebook, which only ever tested
    `== 'E'` / `== 'V'` (never `!= 'E'`). Genes with essentiality ==
    `Not_determined` (no deletion_library_categories.xlsx call for that gene)
    are EXCLUDED from both buckets here (previously an earlier draft folded
    them into "non_essential" via `!= 'E'`, which silently diverged from the
    `_HIST_ROW_QUERIES` == 'V' filter used by plot_dr_dl_histograms and
    produced inconsistent totals between coverage_stats.tsv and the PDF).
    """
    essential = gene_result[gene_result["essentiality"] == "E"]
    non_essential = gene_result[gene_result["essentiality"] == "V"]
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


def _compute_category_coverage(gene_result: pd.DataFrame, column: str) -> dict[str, dict[str, int]]:
    """Split compute_gene_coverage by every non-null value of `column`.

    Returns a dict mapping each value to its coverage stats
    (total/covered/not_covered). Shared by the per-column category breakdowns
    (characterisation_status, deletion_viability, essentiality) that all feed
    build_stats_table's category rows.
    """
    if column not in gene_result.columns:
        logger.warning(f"{column} column not found in gene_result")
        return {}

    result = {}
    value_counts = gene_result[column].value_counts()
    logger.info(f"Computing coverage for {len(value_counts)} {column} categories")

    for value in value_counts.index:
        if pd.isna(value):
            continue
        subset = gene_result[gene_result[column] == value]
        result[value] = compute_gene_coverage(subset)

    return result


def compute_characterisation_status_coverage(gene_result: pd.DataFrame) -> dict[str, dict[str, int]]:
    """Split compute_gene_coverage by characterisation_status values.

    Returns a dict mapping each characterisation_status value to its coverage
    stats (total/covered/not_covered). Only includes protein-coding genes that
    have a non-null characterisation_status annotation.
    """
    return _compute_category_coverage(gene_result, "characterisation_status")


def compute_deletion_viability_coverage(gene_result: pd.DataFrame) -> dict[str, dict[str, int]]:
    """Split compute_gene_coverage by deletion_viability values.

    deletion_viability has 4 categories (viable/inviable/depends_on_conditions/
    unknown), all sourced from gene_metadata for the full protein-coding gene
    universe, so every category is represented here (no nulls to skip).
    """
    return _compute_category_coverage(gene_result, "deletion_viability")


def compute_essentiality_category_coverage(gene_result: pd.DataFrame) -> dict[str, dict[str, int]]:
    """Split compute_gene_coverage by essentiality values (E / V / Not_determined).

    Unlike compute_essentiality_coverage (which keeps the E/V two-bucket split
    used by the donut plots and DR/DL histograms, excluding Not_determined
    genes entirely), this gives every essentiality value — including
    Not_determined — its own coverage_stats.tsv row.
    """
    return _compute_category_coverage(gene_result, "essentiality")


def build_detailed_gene_table(gene_result: pd.DataFrame, gene_metadata: pd.DataFrame) -> pd.DataFrame:
    """Build a detailed gene-level table with DIT-HAP data + metadata for all protein-coding genes.

    Returns a table with columns:
    - Systematic ID, Name, product (from metadata)
    - characterisation_status, deletion_viability (from metadata)
    - DR, DL, essentiality (from gene_result, DR/DL NaN if not covered;
      essentiality is never null — "Not_determined" when no deletion-library call exists)
    - coverage_status: "covered" if DR is not NaN, "not_covered" otherwise

    Sorted by characterisation_status (descending by gene count), then by coverage_status, then by DR (desc).
    """
    # Start with full protein-coding gene universe from metadata
    protein_genes = gene_metadata[gene_metadata["feature_type"] == "protein"].copy()

    # Select commonly used metadata columns
    meta_cols = ["systematic_id", "name", "product", "characterisation_status", "deletion_viability"]
    available_meta_cols = [c for c in meta_cols if c in protein_genes.columns]
    base_table = protein_genes[available_meta_cols].copy()
    base_table = base_table.rename(columns={"systematic_id": "Systematic ID", "name": "Name"})

    # Merge with gene_result (left join so uncovered genes remain)
    dit_hap_cols = ["Systematic ID", "DR", "DL", "essentiality"]
    available_dit_hap_cols = [c for c in dit_hap_cols if c in gene_result.columns]

    detailed_table = base_table.merge(
        gene_result[available_dit_hap_cols],
        on="Systematic ID",
        how="left"
    )

    # Add coverage status
    detailed_table["coverage_status"] = detailed_table["DR"].notna().map({True: "covered", False: "not_covered"})

    # Sort: by characterisation_status frequency (most common first), then coverage, then DR
    if "characterisation_status" in detailed_table.columns:
        status_order = detailed_table["characterisation_status"].value_counts().index.tolist()
        detailed_table["_status_rank"] = detailed_table["characterisation_status"].map(
            {s: i for i, s in enumerate(status_order)}
        )
        detailed_table = detailed_table.sort_values(
            ["_status_rank", "coverage_status", "DR"],
            ascending=[True, True, False],
            na_position="last"
        ).drop(columns=["_status_rank"])
    else:
        detailed_table = detailed_table.sort_values(
            ["coverage_status", "DR"],
            ascending=[True, False],
            na_position="last"
        )

    return detailed_table.reset_index(drop=True)


def write_detailed_gene_excel(detailed_table: pd.DataFrame, output_path: Path) -> None:
    """Write detailed gene table to Excel with multiple sheets: one for all genes, then one per characterisation_status.

    Sheets:
    - "All genes": complete table (5,126 genes)
    - "biological role published", "biological role inferred", etc.: one sheet per characterisation_status category
    """
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        # Sheet 1: All genes
        detailed_table.to_excel(writer, sheet_name="All genes", index=False)

        # Sheets 2+: One per characterisation_status category
        if "characterisation_status" in detailed_table.columns:
            status_counts = detailed_table["characterisation_status"].value_counts()
            for status in status_counts.index:
                if pd.isna(status):
                    continue
                subset = detailed_table[detailed_table["characterisation_status"] == status]
                # Excel sheet names are limited to 31 characters
                sheet_name = str(status)[:31]
                subset.to_excel(writer, sheet_name=sheet_name, index=False)

    logger.info(f"Wrote detailed gene Excel: {len(detailed_table):,} genes, {len(status_counts)} characterisation_status sheets")


# =============================================================================
# STATS TABLE ASSEMBLY
# =============================================================================
def build_stats_table(
    insertion_coverage: dict[str, int],
    gene_coverage: dict[str, int],
    essentiality_coverage: dict[str, dict[str, int]],
    per_chromosome: pd.DataFrame,
    characterisation_status_coverage: dict[str, dict[str, int]] | None = None,
    deletion_viability_coverage: dict[str, dict[str, int]] | None = None,
    essentiality_category_coverage: dict[str, dict[str, int]] | None = None,
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

    # Per-category coverage rows: characterisation_status (arbitrary # of
    # categories), deletion_viability (4: viable/inviable/depends_on_conditions/
    # unknown), essentiality (3: E/V/Not_determined — the full breakdown, unlike
    # essentiality_coverage's essential/non_essential rows above which exclude
    # Not_determined). Each dict contributes one row per category, prefixed so
    # coverage_dicts_from_stats_table can recover which breakdown a row belongs to.
    for prefix, coverage in (
        ("characterisation_", characterisation_status_coverage),
        ("deletion_viability_", deletion_viability_coverage),
        ("essentiality_", essentiality_category_coverage),
    ):
        if not coverage:
            continue
        for category, counts in coverage.items():
            rows.append({
                "metric": "gene",
                "category": f"{prefix}{category}",
                "total": counts["total"],
                "covered": counts["covered"],
                "not_covered": counts["not_covered"],
            })

    stats = pd.DataFrame(rows)

    # Percent columns (of total), rounded to 1 decimal. total == 0 -> NaN rather
    # than a divide-by-zero (no category should be empty, but stay defensive).
    total = stats["total"].replace(0, np.nan)
    stats["covered_pct"] = (stats["covered"] / total * 100).round(1)
    stats["not_covered_pct"] = (stats["not_covered"] / total * 100).round(1)

    # Group rows by metric (in first-appearance order: insertion, then gene), then
    # sort by total descending within each group. Stable mergesort keeps the metric
    # grouping intact while ordering each group's rows biggest-first.
    metric_rank = {m: i for i, m in enumerate(stats["metric"].drop_duplicates())}
    stats = (
        stats.assign(_metric_rank=stats["metric"].map(metric_rank))
        .sort_values(["_metric_rank", "total"], ascending=[True, False], kind="stable")
        .drop(columns="_metric_rank")
        .reset_index(drop=True)
    )

    return stats


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


def _grid_shape(n: int, ncols: int) -> tuple[int, int]:
    """Rows/cols for an n-panel grid at a fixed column count (ceil-divide for the last partial row)."""
    nrows = (n + ncols - 1) // ncols
    return nrows, ncols


def plot_characterisation_status_donuts(
    characterisation_status_coverage: dict[str, dict[str, int]],
    ncols: int = 4,
) -> plt.Figure:
    """One covered/not-covered donut per characterisation_status category, in a grid.

    Categories are drawn in the dict's own order (compute_characterisation_status_coverage
    yields them most-frequent-first). Any trailing empty grid cells are hidden.
    """
    statuses = list(characterisation_status_coverage)
    nrows, ncols = _grid_shape(len(statuses), ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(AX_WIDTH * ncols, AX_HEIGHT * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, status in zip(axes, statuses):
        counts = characterisation_status_coverage[status]
        donut_chart(
            values=[counts["covered"], counts["not_covered"]],
            labels=["Covered", "Not covered"],
            colors=_CHARACTERISATION_COVERAGE_COLORS,
            center_text=f"Total\n{counts['total']:,}\ngenes",
            ax=ax,
        )
        # Wrap long status labels so titles don't collide with neighbours.
        ax.set_title(_wrap_label(status))

    for ax in axes[len(statuses):]:
        ax.axis("off")

    fig.suptitle("Gene coverage by characterisation status", y=1.02)
    fig.tight_layout(h_pad=1, w_pad=1)
    return fig


def plot_characterisation_status_histograms(
    gene_result: pd.DataFrame,
    ncols: int = 2,
) -> plt.Figure:
    """DR + DL distribution histograms per characterisation_status category.

    One row per category, two columns (DR, DL) — mirroring plot_dr_dl_histograms'
    per-feature bins/x-limits. Categories are ordered most-frequent-first.
    """
    if "characterisation_status" not in gene_result.columns:
        logger.warning("characterisation_status column not found; skipping per-category histograms")
        fig, ax = plt.subplots()
        ax.axis("off")
        return fig

    statuses = gene_result["characterisation_status"].value_counts().index.tolist()
    features = ["DR", "DL"]
    nrows = len(statuses)
    fig, axes = plt.subplots(nrows, len(features), figsize=(AX_WIDTH * len(features), AX_HEIGHT * nrows))
    axes = np.atleast_2d(axes)

    for row, status in enumerate(statuses):
        subset = gene_result[gene_result["characterisation_status"] == status]
        for col, feature in enumerate(features):
            ax = axes[row, col]
            bins, xlim = (_DR_BINS, _DR_XLIM) if feature == "DR" else (_DL_BINS, _DL_XLIM)
            data = subset[feature].dropna()
            ax.hist(data, bins=bins, rwidth=0.9, color=_CHARACTERISATION_COVERAGE_COLORS[0])
            ax.set_xlim(xlim)
            if col == 0:
                ax.set_ylabel(f"{_wrap_label(status)}\nNumber of genes")
            if row == 0:
                ax.set_title(f"{feature} distribution")

    fig.tight_layout()
    return fig


def _wrap_label(label: str, width: int = 20) -> str:
    """Soft-wrap a long category label onto multiple lines for plot titles/axis labels."""
    import textwrap

    return "\n".join(textwrap.wrap(str(label), width=width)) or str(label)
