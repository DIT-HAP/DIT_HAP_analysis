#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Deletion Library Phenotype Verification
==========================================

Per-dataset: merges gene-level DIT-HAP results with the Hayles-2013-derived
deletion library phenotype categories and the curated essentiality
verification table, then reports category-level counts and a DR-vs-category
scatter. Ported from
DIT_HAP_pipeline/workflow/notebooks/compare_with_deletion_library.ipynb.

NOTE: the curated deletion_library_categories.xlsx schema changed after the
source notebook was written — `Updated_Systematic_ID` no longer exists,
`Systematic ID` now holds the current ID directly, and `Category` values
changed (e.g. `WT` -> `WT-like`). merge_deletion_library() accepts either
column name so it works against both the old-schema test fixture and the
current real file (see deletion-library-categories-schema-update memory note).
This is a simplified port: the notebook's altair area-fraction-vs-DR charts,
boxplot/violin comparisons, and per-gene depletion-curve PDFs are out of scope
here — only the category breakdown (donut) + DR scatter the task calls for.

Input
-----
- Gene-level fitting_results.tsv (Systematic ID, DR, DL,
  DeletionLibrary_essentiality already native columns — same as coverage.smk).
- resources/curated/deletion_library_categories.xlsx (Systematic ID or
  Updated_Systematic_ID + Category).
- resources/curated/essentiality_verification.csv (systematic_id,
  verification_phenotype, verification_essentiality, ...).

Output
------
- verification_stats.tsv: category counts + verification match/mismatch counts.
- deletion_library_comparison.pdf: donut chart of phenotype categories + DR
  scatter by category.

Usage
-----
    python compare_deletion_library.py \\
        --fitting-results .../release/gene_level/fitting_results.tsv \\
        --deletion-library resources/curated/deletion_library_categories.xlsx \\
        --essentiality-verification resources/curated/essentiality_verification.csv \\
        --output-stats results/verification/{dataset}/verification_stats.tsv \\
        --output-figures results/verification/{dataset}/deletion_library_comparison.pdf

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
from workflow.src.plotting.style import (  # noqa: E402
    AX_HEIGHT,
    AX_WIDTH,
    CATEGORY_COLOR_MAP,
    DONUT_COLOR_MAP,
)


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Byte-faithful to the source notebook's simplified_verification_result: some
# curated verification_phenotype rows carry a compound label (multiple
# phenotypes observed) or a growth-condition caveat; all three collapse to
# the plain "E" essentiality-verified bucket.
_VERIFICATION_PHENOTYPE_SIMPLIFY = {
    "E,small colonies": "E",
    "E,WT": "E",
    "Leu-condition": "E",
}

# Category display order for the donut chart, byte-faithful to the notebook's
# label_orders (categories not present in the data are silently skipped).
_CATEGORY_ORDER = ["spores", "germinated", "microcolonies", "E", "very small colonies", "small colonies", "WT"]

# Legacy -> current metric column names, same quirk as coverage.smk's
# compute_coverage_stats.load_gene_level / clustering's candidates.load_and_annotate:
# some releases' gene-level fitting_results.tsv still ship the pre-rename um/lam
# headers instead of DR/DL.
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}

# Display-only alias for the current deletion_library_categories.xlsx schema:
# CATEGORY_COLOR_MAP / DONUT_COLOR_MAP only know the notebook's original,
# single-phenotype vocabulary, but the current curated file also has (a) the
# straight rename "WT" -> "WT-like", and (b) compound multi-phenotype labels
# that didn't exist when the notebook was written (e.g. a gene scored as both
# "spores" and "germinated" for the same category). Each compound label is
# folded into the single existing bucket judged most representative of it, so
# every gene keeps a legible color/x-position instead of falling back to gray
# or getting dropped — see the per-key comments below for the specific choice.
# compute_category_stats / merge_deletion_library / apply_category_with_essentiality
# are untouched by this and keep reporting the raw Category string; only
# plot_category_donut / plot_dr_scatter_by_category consult this mapping.
_CATEGORY_DISPLAY_ALIASES = {
    "WT-like": "WT",
    # First-listed phenotype is the one that determines the merged bucket for
    # all four compound "spores, ..." labels below - "spores" is the most
    # severe/earliest-observed phenotype PomBase records for these genes, so
    # it's the most informative single label to keep.
    "spores, germinated": "spores",
    "spores, germinated, divided or microcolonies": "spores",
    "spores, miscellaneous": "spores",
    # "germinated, divided or microcolonies" -> "germinated": same logic,
    # germination is the earliest-observed phenotype in this compound label.
    "germinated, divided or microcolonies": "germinated",
    # "microcolonies, small colonies" -> "small colonies": the second term is
    # the more specific/severe of the two (a true "small colony" is a finer
    # distinction than the broader "microcolonies" bucket).
    "microcolonies, small colonies": "small colonies",
}


def _display_category(category: str) -> str:
    """Map a raw Category value to the color-map key used for plotting (see _CATEGORY_DISPLAY_ALIASES)."""
    return _CATEGORY_DISPLAY_ALIASES.get(category, category)


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class VerificationConfig:
    """Inputs and outputs for the deletion library verification analysis."""
    fitting_results: Path
    deletion_library: Path
    essentiality_verification: Path
    output_stats: Path
    output_figures: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.fitting_results, self.deletion_library, self.essentiality_verification]:
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


def load_deletion_library(deletion_library_path: Path) -> pd.DataFrame:
    """Load the curated deletion library xlsx, keeping just the ID + Category columns.

    Handles both schemas: the old file's `Updated_Systematic_ID` key (also used
    by the unit test fixture) and the current file's `Systematic ID` key (see
    module docstring / deletion-library-categories-schema-update memory note).
    """
    deletion_library = pd.read_excel(deletion_library_path)
    id_col = "Updated_Systematic_ID" if "Updated_Systematic_ID" in deletion_library.columns else "Systematic ID"
    return deletion_library[[id_col, "Category"]]


def load_essentiality_verification(essentiality_verification_path: Path) -> pd.DataFrame:
    """Load + rename the curated essentiality verification table to the notebook's column names.

    Simplifies compound `verification_phenotype` values (multi-phenotype or
    growth-condition-caveated) down to plain "E" (see
    _VERIFICATION_PHENOTYPE_SIMPLIFY), then drops rows missing either the
    phenotype or essentiality call.
    """
    verification = pd.read_csv(essentiality_verification_path).rename(
        columns={
            "systematic_id": "Systematic ID",
            "verification_phenotype": "Verification result",
            "verification_essentiality": "Verified essentiality",
        }
    )[["Systematic ID", "Verification result", "Verified essentiality"]]
    verification["Verification result"] = verification["Verification result"].replace(
        _VERIFICATION_PHENOTYPE_SIMPLIFY
    )
    return verification.dropna(subset=["Verification result", "Verified essentiality"])


# =============================================================================
# CORE LOGIC — merge + stats (unit-tested)
# =============================================================================
def merge_deletion_library(gene_result: pd.DataFrame, deletion_library: pd.DataFrame) -> pd.DataFrame:
    """Left-merge gene-level results with deletion library categories on Systematic ID.

    Accepts either the old schema (`Updated_Systematic_ID` key) or the current
    schema (`Systematic ID` key) for `deletion_library` — see module docstring.
    """
    if "Updated_Systematic_ID" in deletion_library.columns:
        return gene_result.merge(
            deletion_library,
            left_on="Systematic ID",
            right_on="Updated_Systematic_ID",
            how="left",
        ).drop(columns=["Updated_Systematic_ID"])
    if "Systematic ID" in deletion_library.columns:
        return gene_result.merge(deletion_library, on="Systematic ID", how="left")
    raise KeyError(
        "deletion_library must contain a 'Systematic ID' or 'Updated_Systematic_ID' column"
    )


def apply_category_with_essentiality(row: pd.Series) -> str:
    """Append an "(E)" suffix to 'small colonies' when DeletionLibrary_essentiality == 'E'.

    Byte-faithful to the notebook's Category_with_essentiality column: flags
    the subset of nominally-"small colonies" genes that the deletion library
    itself also calls essential, without altering any other category label.
    """
    if row["Category"] == "small colonies" and row["DeletionLibrary_essentiality"] == "E":
        return f"{row['Category']} (E)"
    return row["Category"]


def compute_category_stats(merged: pd.DataFrame) -> pd.DataFrame:
    """Count genes per deletion-library Category in the merged frame."""
    return (
        merged.groupby("Category", dropna=False)
        .size()
        .reset_index(name="count")
        .rename(columns={"Category": "category"})
    )


def compute_category_with_essentiality_stats(merged: pd.DataFrame) -> pd.DataFrame:
    """Count genes per Category_with_essentiality (see apply_category_with_essentiality).

    Surfaces the "small colonies (E)" split that apply_category_with_essentiality
    computes but compute_category_stats' plain Category grouping doesn't show
    on its own — otherwise Category_with_essentiality would be computed and
    never read anywhere.
    """
    return (
        merged.groupby("Category_with_essentiality", dropna=False)
        .size()
        .reset_index(name="count")
        .rename(columns={"Category_with_essentiality": "category"})
    )


def merge_essentiality_verification(merged: pd.DataFrame, essentiality_verification: pd.DataFrame) -> pd.DataFrame:
    """Left-merge in the curated Verification result / Verified essentiality columns."""
    return merged.merge(essentiality_verification, on="Systematic ID", how="left")


def compute_verification_match_stats(merged_with_verification: pd.DataFrame) -> dict[str, int]:
    """Compare curated 'Verified essentiality' against 'DeletionLibrary_essentiality' where both are known."""
    known = merged_with_verification.dropna(subset=["Verified essentiality", "DeletionLibrary_essentiality"])
    match = int((known["Verified essentiality"] == known["DeletionLibrary_essentiality"]).sum())
    return {
        "verified_total": len(known),
        "match": match,
        "mismatch": len(known) - match,
    }


# =============================================================================
# STATS TABLE ASSEMBLY
# =============================================================================
def build_stats_table(
    category_stats: pd.DataFrame,
    category_with_essentiality_stats: pd.DataFrame,
    verification_stats: dict[str, int],
) -> pd.DataFrame:
    """Flatten category counts + verification match/mismatch counts into one long-form table."""
    rows = [
        {"metric": "category_count", "category": row["category"], "count": row["count"]}
        for _, row in category_stats.iterrows()
    ]
    rows += [
        {"metric": "category_with_essentiality_count", "category": row["category"], "count": row["count"]}
        for _, row in category_with_essentiality_stats.iterrows()
    ]
    for key, count in verification_stats.items():
        rows.append({"metric": "verification", "category": key, "count": count})
    return pd.DataFrame(rows)


# =============================================================================
# PLOTTING
# =============================================================================
def plot_category_donut(category_stats: pd.DataFrame) -> plt.Figure:
    """Donut chart of gene counts per deletion-library phenotype category."""
    raw_categories = set(category_stats["category"])
    ordered = [c for c in raw_categories if _display_category(c) in _CATEGORY_ORDER]
    ordered.sort(key=lambda c: _CATEGORY_ORDER.index(_display_category(c)))
    remaining = [c for c in category_stats["category"] if c not in ordered]
    labels = ordered + remaining
    counts_by_label = category_stats.set_index("category")["count"]
    values = [int(counts_by_label[label]) for label in labels]
    colors = [DONUT_COLOR_MAP.get(_display_category(label), "gray") for label in labels]

    # Anything still falling back to gray is either a genuine merge-failure
    # NaN category or a Category string _CATEGORY_DISPLAY_ALIASES doesn't
    # know about yet (future schema drift) — warn with counts so this isn't
    # silently indistinguishable from a real "gray" bucket in the figure.
    unmapped = [(label, int(counts_by_label[label])) for label, color in zip(labels, colors) if color == "gray"]
    if unmapped:
        logger.warning(
            f"{sum(n for _, n in unmapped):,} genes plotted as gray (unmapped category): "
            + ", ".join(f"{label!r} (n={n})" for label, n in unmapped)
        )

    fig, ax = plt.subplots(figsize=(AX_WIDTH, AX_HEIGHT))
    donut_chart(
        values=values,
        labels=labels,
        colors=colors,
        center_text=f"Total\n{sum(values):,}\ngenes",
        ax=ax,
    )
    ax.set_title("Deletion library phenotype categories")
    fig.tight_layout()
    return fig


def plot_dr_scatter_by_category(merged: pd.DataFrame) -> plt.Figure:
    """Scatter of DR per gene, grouped by category (x-jittered for visibility).

    Per-category n= counts are shown directly on the x-tick labels rather
    than a legend: with _CATEGORY_ORDER plus every distinct compound-label
    alias this can be 10+ categories, and a legend key of that size clutters
    a scatter more than it helps.
    """
    all_categories = merged["Category"].dropna().unique()
    categories = [c for c in all_categories if _display_category(c) in CATEGORY_COLOR_MAP]
    categories.sort(key=lambda c: _CATEGORY_ORDER.index(_display_category(c)))

    excluded = [c for c in all_categories if c not in categories]
    if excluded:
        excluded_genes = merged["Category"].isin(excluded).sum()
        logger.warning(
            f"{excluded_genes:,} genes excluded from DR scatter (category has no display color mapping): "
            + ", ".join(str(c) for c in excluded)
        )

    fig, ax = plt.subplots(figsize=(AX_WIDTH, AX_HEIGHT))
    rng = np.random.default_rng(42)
    tick_labels = []
    for i, category in enumerate(categories):
        dr_values = merged.query("Category == @category")["DR"].dropna()
        jitter = rng.uniform(-0.15, 0.15, size=len(dr_values))
        ax.scatter(
            i + jitter, dr_values,
            alpha=0.5, s=10,
            color=CATEGORY_COLOR_MAP[_display_category(category)],
        )
        tick_labels.append(f"{category}\n(n={len(dr_values)})")
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(tick_labels, rotation=30, ha="right")
    ax.set_ylabel("Depletion Rate (DR)")
    ax.set_title("DR by deletion library category")
    fig.tight_layout()
    return fig


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: VerificationConfig) -> None:
    """Load -> merge -> compute stats -> save TSV + figures."""
    config.validate()

    gene_result = load_gene_level(config.fitting_results)
    deletion_library = load_deletion_library(config.deletion_library)
    essentiality_verification = load_essentiality_verification(config.essentiality_verification)

    merged = merge_deletion_library(gene_result, deletion_library)
    merged["Category_with_essentiality"] = merged.apply(apply_category_with_essentiality, axis=1)
    merged_with_verification = merge_essentiality_verification(merged, essentiality_verification)

    category_stats = compute_category_stats(merged)
    category_with_essentiality_stats = compute_category_with_essentiality_stats(merged)
    verification_stats = compute_verification_match_stats(merged_with_verification)

    stats_table = build_stats_table(category_stats, category_with_essentiality_stats, verification_stats)
    stats_table.to_csv(config.output_stats, sep="\t", index=False)

    fig_donut = plot_category_donut(category_stats)
    fig_scatter = plot_dr_scatter_by_category(merged)

    with PdfPages(config.output_figures) as pdf:
        pdf.savefig(fig_donut, dpi=300, bbox_inches="tight")
        pdf.savefig(fig_scatter, dpi=300, bbox_inches="tight")
    plt.close(fig_donut)
    plt.close(fig_scatter)

    logger.success(
        f"Verification: {len(merged):,} genes across {len(category_stats):,} categories, "
        f"{verification_stats['match']:,}/{verification_stats['verified_total']:,} curated verifications match"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Compare gene-level results against the deletion library phenotype categories")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Gene-level fitting_results.tsv")
    parser.add_argument("--deletion-library", type=Path, required=True, help="Curated deletion_library_categories.xlsx")
    parser.add_argument("--essentiality-verification", type=Path, required=True, help="Curated essentiality_verification.csv")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output verification stats TSV")
    parser.add_argument("--output-figures", type=Path, required=True, help="Output verification figures PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = VerificationConfig(
            fitting_results=args.fitting_results,
            deletion_library=args.deletion_library,
            essentiality_verification=args.essentiality_verification,
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
