"""
Gene Insertion Coverage — Core Logic
=====================================

Shared constants, loaders, coverage computations, stats-table assembly, and the
render-ready frame builders for the coverage figures. Ported from
DIT_HAP_pipeline/workflow/notebooks/gene_coverage_analysis.ipynb and factored
out of the original single-script port so the stage can be split into
independent Snakemake rules (prepare -> compute stats / plot figures), each
re-runnable on its own.

Drawing itself belongs to workflow/src/figure_render/; the builders at the
bottom of this module only reshape numbers into the frames those renderers take.

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
  prepare_coverage_data.py: the gene universe and its annotation columns
  (characterisation_status, FYPOviability, deletion_essentiality) come from the
  gene annotation reference instead, so coverage reads the same values — under
  the same column names — as the rest of the analysis (see
  prepare_coverage_data.py's run()).

Usage
-----
    from coverage.core import (
        load_insertion_level, resolve_duplicate_annotations,
        compute_insertion_coverage, compute_gene_coverage,
        compute_per_chromosome_insertion_coverage,
        compute_category_coverage,
        build_stats_table, coverage_dicts_from_stats_table,
        composition_frame, dimension_coverage_frame, insertion_placement_frame,
        dr_dl_histogram_frame, DIMENSION_LABELS,
    )
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from collections.abc import Mapping
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
from release_schema import IN_GENE_FILTER, read_gene_level


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# DR/DL histogram bin edges, byte-faithful to the notebook's "DR DL Histogram"
# cell. Passed to the histogram renderer as explicit edges, so they also set
# each panel's x range (the notebook's own xlim).
#
# Mirrored 2026-09-17 when upstream flipped the DR sign convention (negative =
# depleted). The DR edges are the exact negation of the old ones, so the axis
# still truncates the extreme-depletion tail rather than the WT shoulder — the
# old [-0.2, 1.45] cut the top 0.1%, the new [-1.45, 0.2] cuts the bottom ~1%.
# Leaving them un-mirrored would have silently dropped 38% of genes (all the
# depleted ones) out of the histogram. DL did not flip.
DR_BINS = np.arange(-1.45, 0.25, 0.05)
DL_BINS = np.arange(0, 15, 0.5)

# The three per-gene annotation dimensions every coverage breakdown is computed
# for: (gene_result column, display labels for that column's values). The column
# names are the annotation reference's own (1c_annotate.smk), so a coverage table
# reads the same as the reference and the annotated workbook, with no renaming
# in between. Labels are display-only; keys absent from a map fall back to the
# raw value, so a new category still renders rather than vanishing.
#
# One ordered registry, consumed by everything downstream: compute_coverage_stats
# computes a breakdown per column, build_stats_table prefixes its rows with the
# column name, coverage_dicts_from_stats_table reads them back, and
# plot_coverage_figures renders the SAME two figures per column (coverage
# composition + DR/DL distributions) — so no dimension is analysed differently
# from the others, and adding one is a line here.
DELETION_VIABILITY_LABELS = {
    "viable": "Viable",
    "inviable": "Inviable",
    "depends_on_conditions": "Depends on conditions",
    "unknown": "Unknown",
}
DELETION_ESSENTIALITY_LABELS = {
    "E": "Essential",
    "V": "Non-essential",
    "Not_determined": "Not determined",
}
DIMENSION_LABELS = {
    "characterisation_status": {},
    "FYPOviability": DELETION_VIABILITY_LABELS,
    "deletion_essentiality": DELETION_ESSENTIALITY_LABELS,
}

# Insertion placement is plotted for the three main chromosomes only. The
# telomeric gap (220 insertions), mitochondrial (25) and mating-type region (3)
# are too small for their percentage to mean anything beside chr_I's 42,252 —
# 2 of 3 insertions reads as 66.7%.
PLOTTED_CHROMOSOMES = ("I", "II", "III")

# Donut labels. Gene-level figures show a covered fraction against the
# remainder; the insertion-placement figure shows where insertions landed.
COVERED_LABEL = "Covered"
NOT_COVERED_LABEL = "Not covered"
IN_GENE_LABEL = "In genes"
INTERGENIC_LABEL = "Intergenic"


# =============================================================================
# LOADERS
# =============================================================================
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
def coverage_dicts_from_stats_table(
    stats: pd.DataFrame,
) -> tuple[
    dict[str, int],
    pd.DataFrame,
    dict[str, dict[str, dict[str, int]]],
]:
    """Reconstruct the coverage dicts + per-chromosome table from a coverage_stats.tsv frame.

    Inverse of build_stats_table: lets plot_coverage_figures render donuts from the
    exact numbers compute_coverage_stats wrote, instead of recomputing them from the
    gene_result parquet (which risks figure/table drift if the two paths ever diverge).
    Returns (gene_coverage, per_chromosome, dimension_coverage), the three things the
    figures need; dimension_coverage maps each column in DIMENSION_LABELS to its
    {value: counts} breakdown, recovered from the `<column>_` row prefix.
    """
    def _row(metric: str, category: str) -> pd.Series:
        hit = stats[(stats["metric"] == metric) & (stats["category"] == category)]
        if hit.empty:
            raise ValueError(f"coverage_stats table missing required row: metric={metric!r}, category={category!r}")
        return hit.iloc[0]

    gene = _row("gene", "all")
    gene_coverage = {"total": int(gene["total"]), "covered": int(gene["covered"]), "not_covered": int(gene["not_covered"])}

    # Per-chromosome insertion rows: any insertion row that isn't the "all" summary.
    per_chr_rows = stats[(stats["metric"] == "insertion") & (stats["category"] != "all")]
    per_chromosome = pd.DataFrame({
        "Chr": per_chr_rows["category"].str.replace(r"^chr_", "", regex=True),
        "total": per_chr_rows["total"].astype(int),
        "in_gene": per_chr_rows["covered"].astype(int),
        "intergenic": per_chr_rows["not_covered"].astype(int),
    }).reset_index(drop=True)

    def _dimension_coverage(column: str) -> dict[str, dict[str, int]]:
        prefix = f"{column}_"
        result: dict[str, dict[str, int]] = {}
        rows = stats[(stats["metric"] == "gene") & (stats["category"].str.startswith(prefix))]
        for _, r in rows.iterrows():
            value = r["category"][len(prefix):]
            result[value] = {"total": int(r["total"]), "covered": int(r["covered"]), "not_covered": int(r["not_covered"])}
        return result

    dimension_coverage = {column: _dimension_coverage(column) for column in DIMENSION_LABELS}

    return (
        gene_coverage,
        per_chromosome,
        dimension_coverage,
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


def compute_per_chromosome_insertion_coverage(annotation: pd.DataFrame) -> pd.DataFrame:
    """Per-chromosome in-gene/intergenic insertion counts (Chr is the 1st index level)."""
    rows = []
    for chrom, group in annotation.groupby(level="Chr"):
        counts = compute_insertion_coverage(group)
        rows.append({"Chr": chrom, **counts})
    return pd.DataFrame(rows).sort_values("Chr").reset_index(drop=True)


def compute_category_coverage(gene_result: pd.DataFrame, column: str) -> dict[str, dict[str, int]]:
    """Split compute_gene_coverage by every non-null value of `column`.

    Returns a dict mapping each value to its coverage stats
    (total/covered/not_covered), most common value first. Called once per column
    in DIMENSION_LABELS, so all three annotation dimensions get an identical
    breakdown; `column` missing from the table is a warning, not a crash, so a
    reference built without a block still produces the other figures.
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


def build_detailed_gene_table(gene_result: pd.DataFrame) -> pd.DataFrame:
    """Build a detailed gene-level table with DIT-HAP data + annotation for all protein-coding genes.

    `gene_result` is already the full protein-coding gene universe, carrying its annotation
    columns (Name, product, characterisation_status, FYPOviability, deletion_essentiality) from
    the gene annotation reference — see prepare_coverage_data.

    Returns a table with columns:
    - Systematic ID, Name, product, characterisation_status, FYPOviability
    - DR, DL (NaN if not covered)
    - deletion_essentiality (never null — "Not_determined" when no deletion-library call exists)
    - coverage_status: "covered" if DR is not NaN, "not_covered" otherwise

    Sorted by characterisation_status (descending by gene count), then by coverage_status,
    then by DR ascending — most depleted first, since negative DR now means depleted.
    """
    detail_cols = [
        "Systematic ID", "Name", "product", "characterisation_status",
        "FYPOviability", "DR", "DL", "deletion_essentiality",
    ]
    detailed_table = gene_result[[c for c in detail_cols if c in gene_result.columns]].copy()

    # Add coverage status
    detailed_table["coverage_status"] = detailed_table["DR"].notna().map({True: "covered", False: "not_covered"})

    # Sort: by characterisation_status frequency (most common first), then coverage, then DR
    # ascending (most depleted first — negative DR is the depleted end).
    if "characterisation_status" in detailed_table.columns:
        status_order = detailed_table["characterisation_status"].value_counts().index.tolist()
        detailed_table["_status_rank"] = detailed_table["characterisation_status"].map(
            {s: i for i, s in enumerate(status_order)}
        )
        detailed_table = detailed_table.sort_values(
            ["_status_rank", "coverage_status", "DR"],
            ascending=[True, True, True],
            na_position="last"
        ).drop(columns=["_status_rank"])
    else:
        detailed_table = detailed_table.sort_values(
            ["coverage_status", "DR"],
            ascending=[True, True],
            na_position="last"
        )

    return detailed_table.reset_index(drop=True)


# deletion_essentiality's raw E/V/Not_determined values aren't self-descriptive as sheet
# tabs, so map them to the same essential/non_essential naming already used by
# essentiality_coverage's bucket names (FYPOviability's values are used as-is).
_ESSENTIALITY_SHEET_NAMES = {"E": "essential", "V": "non_essential", "Not_determined": "essentiality_not_determined"}


def write_detailed_gene_excel(detailed_table: pd.DataFrame, output_path: Path) -> None:
    """Write detailed gene table to Excel with multiple sheets.

    Sheets:
    - "All genes": complete table (5,126 genes)
    - one sheet per characterisation_status category (e.g. "biological role published")
    - one sheet per deletion_essentiality category ("essential" / "non_essential" / "essentiality_not_determined")
    - one sheet per FYPOviability category ("viable" / "inviable" / "depends_on_conditions" / "unknown")
    """
    sheet_count = 1
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        # Sheet 1: All genes
        detailed_table.to_excel(writer, sheet_name="All genes", index=False)

        # One sheet per characterisation_status category
        if "characterisation_status" in detailed_table.columns:
            for status in detailed_table["characterisation_status"].value_counts().index:
                if pd.isna(status):
                    continue
                subset = detailed_table[detailed_table["characterisation_status"] == status]
                # Excel sheet names are limited to 31 characters
                subset.to_excel(writer, sheet_name=str(status)[:31], index=False)
                sheet_count += 1

        # One sheet per deletion_essentiality category
        if "deletion_essentiality" in detailed_table.columns:
            for value, sheet_name in _ESSENTIALITY_SHEET_NAMES.items():
                subset = detailed_table[detailed_table["deletion_essentiality"] == value]
                if subset.empty:
                    continue
                subset.to_excel(writer, sheet_name=sheet_name[:31], index=False)
                sheet_count += 1

        # One sheet per FYPOviability category
        if "FYPOviability" in detailed_table.columns:
            for viability in detailed_table["FYPOviability"].value_counts().index:
                if pd.isna(viability):
                    continue
                subset = detailed_table[detailed_table["FYPOviability"] == viability]
                subset.to_excel(writer, sheet_name=str(viability)[:31], index=False)
                sheet_count += 1

    logger.info(f"Wrote detailed gene Excel: {len(detailed_table):,} genes, {sheet_count} sheets")


# =============================================================================
# STATS TABLE ASSEMBLY
# =============================================================================
def build_stats_table(
    insertion_coverage: dict[str, int],
    gene_coverage: dict[str, int],
    per_chromosome: pd.DataFrame,
    dimension_coverage: Mapping[str, Mapping[str, Mapping[str, int]]] | None = None,
) -> pd.DataFrame:
    """Flatten all coverage dicts into one long-form stats table.

    `dimension_coverage` maps each column in DIMENSION_LABELS to its
    {value: counts} breakdown; every row is labeled `<column>_<value>`, which is
    both self-describing in the TSV and exactly the key
    coverage_dicts_from_stats_table strips back off. The E/V split lives in there
    as `deletion_essentiality_E` / `_V`, so it is not repeated as its own rows.
    """
    rows = [
        {"metric": "insertion", "category": "all", "total": insertion_coverage["total"],
         "covered": insertion_coverage["in_gene"], "not_covered": insertion_coverage["intergenic"]},
        {"metric": "gene", "category": "all", "total": gene_coverage["total"],
         "covered": gene_coverage["covered"], "not_covered": gene_coverage["not_covered"]},
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

    # Per-category coverage rows, one per value of each annotation dimension
    # (characterisation_status's categories, FYPOviability's 4, deletion_essentiality's
    # 3 — the full breakdown, Not_determined included).
    for column, coverage in (dimension_coverage or {}).items():
        for category, counts in coverage.items():
            rows.append({
                "metric": "gene",
                "category": f"{column}_{category}",
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
# RENDER-READY FRAMES
# =============================================================================
# Drawing belongs to workflow/src/figure_render/: every coverage breakdown goes
# through render_composition_figure (percentage bar + one part/whole donut per
# category), and the DR/DL distributions through render_grouped_histogram_figure.
# These builders only reshape numbers into the frames those renderers take —
# every value still comes from coverage_dicts_from_stats_table, so a figure can
# never disagree with the stats table it was drawn from.

def composition_frame(
    category_coverage: Mapping[str, Mapping[str, int]],
    *,
    labels: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Reshape {category: {covered, not_covered, total}} into a composition frame.

    ``labels`` renames categories for display; unlisted keys pass through as
    they are. ``covered_pct`` is rounded to 1 decimal, matching
    build_stats_table, so the drawn bar labels read identically to
    coverage_stats.tsv.
    """
    labels = labels or {}
    rows = []
    for category, counts in category_coverage.items():
        total = int(counts["total"])
        covered = int(counts["covered"])
        rows.append({
            "category": labels.get(category, category),
            "covered": covered,
            "not_covered": int(counts["not_covered"]),
            "covered_pct": round(covered / total * 100, 1) if total else float("nan"),
        })
    return pd.DataFrame(rows)


def dimension_coverage_frame(
    gene_coverage: Mapping[str, int],
    category_coverage: Mapping[str, Mapping[str, int]],
    labels: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Every gene, then one row per category of one annotation dimension.

    "All genes" and the category rows partition the same 5,126 genes, so the bars
    read as one decomposition rather than N unrelated numbers. Built for each
    column in DIMENSION_LABELS, so all three dimensions' composition figures have
    the same shape.
    """
    return composition_frame({"All genes": gene_coverage, **category_coverage}, labels=labels)


def insertion_placement_frame(per_chromosome: pd.DataFrame) -> pd.DataFrame:
    """In-gene vs intergenic insertion share for the three main chromosomes.

    Small regions (telomeric gap, mitochondrial, mating type) are dropped — see
    PLOTTED_CHROMOSOMES.
    """
    rows = per_chromosome[per_chromosome["Chr"].isin(PLOTTED_CHROMOSOMES)]
    if rows.empty:
        raise ValueError(f"None of {PLOTTED_CHROMOSOMES} found in the per-chromosome table")
    coverage = {
        f"chr_{row['Chr']}": {
            "covered": int(row["in_gene"]),
            "not_covered": int(row["intergenic"]),
            "total": int(row["total"]),
        }
        for _, row in rows.iterrows()
    }
    return composition_frame(coverage)


def dr_dl_histogram_frame(
    gene_result: pd.DataFrame,
    feature: str,
    column: str,
    labels: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Long-form frame for one DR/DL histogram figure: one row per (stratum, gene).

    Strata are "All genes" plus one per value of one annotation dimension
    (`column`), in the same most-common-first order compute_category_coverage
    uses — so the panels line up with that dimension's composition figure. A gene
    contributes to every stratum it belongs to, so an essential gene appears under
    both "All genes" and "Essential": that is what makes the panels comparable.
    """
    if feature not in gene_result.columns:
        raise ValueError(f"{feature!r} not in the gene_result table")
    if column not in gene_result.columns:
        raise ValueError(f"{column!r} not in the gene_result table; cannot stratify")

    labels = labels or {}
    values = [v for v in gene_result[column].value_counts().index if not pd.isna(v)]

    parts = [pd.DataFrame({"stratum": "All genes", feature: gene_result.loc[gene_result[column].notna(), feature].to_numpy()})]
    for value in values:
        mask = gene_result[column] == value
        parts.append(
            pd.DataFrame({"stratum": labels.get(value, str(value)), feature: gene_result.loc[mask, feature].to_numpy()})
        )
    return pd.concat(parts, ignore_index=True)
