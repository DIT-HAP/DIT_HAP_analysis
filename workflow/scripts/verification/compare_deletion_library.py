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
The notebook's matplotlib §4-5 plotting is migrated here: boxplot+violin DR
comparisons, four critical-gene outlier groups (each with a boxplot, a
verification-composition donut, and a review TSV), and DIT-HAP-vs-gRNA
depletion curves. Only the altair interactive charts (area-fraction-vs-DR)
remain notebook-only — they produce HTML and are exploratory. The extra
outputs are optional: omitting the --output-boxplots/... flags falls back to
the original donut + DR-scatter two-output invocation.

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
- verification_boxplots.pdf (optional): basic + critical-group boxplot/violin
  and per-group verification donuts.
- verification_depletion_curves.pdf (optional): per critical group, DIT-HAP
  (+gRNA) depletion curves for that group's outlier genes.
- critical_genes_{group}.tsv (optional): per-group gene detail for review.

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
from workflow.src.plotting.generic import boxplot_with_violinplot, donut_chart  # noqa: E402
from workflow.src.plotting.gene_level import (  # noqa: E402
    DIT_HAP_GENERATIONS,
    GRNA_GENERATIONS,
    plot_gene_depletion_curve,
)
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


# Basic-boxplot category selection, byte-faithful to the notebook's
# selected_categories (cell 6). Grouping is by Category_with_essentiality after
# restricting to these canonical categories.
_BASIC_BOXPLOT_CATEGORIES = ["spores", "germinated", "microcolonies", "very small colonies", "small colonies", "WT"]

# The four "critical gene" outlier groups (notebook §4.2-4.4). Each filter runs
# against the canonicalized `cat_canon` column (so the notebook's raw labels
# like 'WT' / 'small colonies' still match the current WT-like / compound
# schema). `sort` orders the outlier gene list by DR: WT->nonWT / small->E look
# at the highest-DR outliers first (desc), E->V at the lowest (asc).
_CRITICAL_GROUPS = {
    "WT2nonWT": {
        "filter": "cat_canon == 'WT' and DR > 0.35",
        "sort": "desc",
    },
    "scE2E": {
        "filter": "cat_canon == 'small colonies' and DR > 0.75 and DeletionLibrary_essentiality == 'E'",
        "sort": "desc",
    },
    "sc2E": {
        "filter": "cat_canon == 'small colonies' and DR > 0.75 and DeletionLibrary_essentiality != 'E'",
        "sort": "desc",
    },
    "E2V": {
        "filter": "cat_canon in ['spores', 'germinated', 'microcolonies'] and DR < 0.35",
        "sort": "asc",
    },
}

# Verification-result bucket order for the critical-group boxplots/donuts,
# byte-faithful to the notebook's prepare_verification_data loop + label_orders.
_VERIFICATION_BUCKET_ORDER = [
    "spores", "germinated", "microcolonies", "E", "E (tiny colonies)",
    "very small colonies", "small colonies", "WT",
]

# gpd1: the source notebook manually appends this gene to the simplified
# verification table (cell 3) as an E/E verified call — it wasn't in the curated
# verification file then and still isn't now, but it IS a real release gene, so
# the append stays meaningful. Kept byte-faithful.
_MANUAL_VERIFICATION_ROWS = pd.DataFrame(
    {"Systematic ID": ["SPBC215.05"], "Verification result": ["E"], "Verified essentiality": ["E"]}
)


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
    # New (notebook §4-5) outputs. gene_timepoints / grna_timepoints feed the
    # depletion curves; grna_timepoints is optional (HD-only, see verification.smk).
    gene_timepoints: Path | None = None
    grna_timepoints: Path | None = None
    output_boxplots: Path | None = None
    output_depletion_curves: Path | None = None
    output_critical_genes_dir: Path | None = None

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.fitting_results, self.deletion_library, self.essentiality_verification]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.gene_timepoints is not None and not self.gene_timepoints.exists():
            raise ValueError(f"Required input not found: {self.gene_timepoints}")
        outputs = [self.output_stats, self.output_figures, self.output_boxplots, self.output_depletion_curves]
        for out in outputs:
            if out is not None:
                out.parent.mkdir(parents=True, exist_ok=True)
        if self.output_critical_genes_dir is not None:
            self.output_critical_genes_dir.mkdir(parents=True, exist_ok=True)


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
    verification = verification.dropna(subset=["Verification result", "Verified essentiality"])
    # Append the notebook's manual gpd1 (SPBC215.05) E/E call (see
    # _MANUAL_VERIFICATION_ROWS). ignore_index keeps a clean RangeIndex.
    return pd.concat([verification, _MANUAL_VERIFICATION_ROWS], ignore_index=True)


def load_essentiality_verification_full(essentiality_verification_path: Path) -> pd.DataFrame:
    """Load the curated verification table keeping ALL columns (area day3-6, comments, ...).

    Feeds build_final_merged / the per-critical-group review TSVs, which need the
    colony-area measurements that load_essentiality_verification drops. The
    `verification_phenotype` values are NOT simplified here (the raw label is
    what a human reviewer wants to see).
    """
    return pd.read_csv(essentiality_verification_path)


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
# CORE LOGIC — critical-gene analysis (unit-tested)
# =============================================================================
def canonicalize_category(merged: pd.DataFrame) -> pd.DataFrame:
    """Add a `cat_canon` column folding schema-drifted Category labels to notebook vocabulary.

    The notebook's outlier filters use the original single-phenotype labels
    ('WT', 'small colonies', 'spores', ...), but the current curated file ships
    'WT-like' and compound labels. `cat_canon` = `_display_category(Category)`
    lets those filters keep matching (see _CRITICAL_GROUPS).
    """
    out = merged.copy()
    out["cat_canon"] = out["Category"].map(lambda c: _display_category(c) if isinstance(c, str) else c)
    return out


def build_final_merged(
    merged: pd.DataFrame,
    verification_full: pd.DataFrame,
) -> pd.DataFrame:
    """Right-join gene-level+category data with the FULL verification table (area columns kept).

    Reconstructs the notebook's `final_merged`: one row per curated-verification
    gene, carrying DR/DL/FYPOviability/DeletionLibrary_essentiality/Category plus
    the raw verification phenotype/essentiality and colony-area day3-6 columns.
    Feeds the per-critical-group review TSVs. Genes verified as essential ('E')
    but missing a day-3 area are zero-filled for the area columns, byte-faithful
    to the notebook (they were confirmed dead, i.e. zero colony area).
    """
    area_cols = [c for c in verification_full.columns if "area" in c]
    final = merged.merge(
        verification_full.rename(columns={"systematic_id": "Systematic ID"}),
        on="Systematic ID",
        how="right",
    )
    e_missing_day3 = final.query(
        "verification_essentiality == 'E' and median_area_day3.isna()",
        engine="python",
    ).index
    final.loc[e_missing_day3, area_cols] = 0
    return final


def prepare_verification_data(
    merged: pd.DataFrame,
    final_merged: pd.DataFrame,
    simplified_verification: pd.DataFrame,
    outlier_filter: str,
    sort: str = "desc",
) -> tuple[dict[str, list[float]], pd.DataFrame]:
    """Bucket a group's outliers by verification result; return {bucket: [DR...]} + gene detail.

    Ported from the notebook's prepare_verification_data. Selects outliers via
    `outlier_filter` (run against `merged`, which must already have `cat_canon`),
    crosses them with the simplified verification table, and buckets into
    {"Not verified": [...], <verified category>: [...]} preserving
    _VERIFICATION_BUCKET_ORDER. Each bucket's value is the list of member DR
    values (for the boxplot); bucket size drives the donut. The second return is
    the per-gene detail frame (from final_merged) tagged with its bucket, for
    the review TSV.
    """
    ascending = sort == "asc"
    outliers = (
        merged.query(outlier_filter, engine="python")
        .sort_values("DR", ascending=ascending)["Systematic ID"]
        .unique()
        .tolist()
    )
    verified = simplified_verification[simplified_verification["Systematic ID"].isin(outliers)]
    verified_genes = set(verified["Systematic ID"])
    missing = [g for g in outliers if g not in verified_genes]

    buckets: dict[str, list[str]] = {}
    if missing:
        buckets["Not verified"] = missing
    for category in _VERIFICATION_BUCKET_ORDER:
        genes = verified.loc[verified["Verification result"] == category, "Systematic ID"].unique().tolist()
        if genes:
            buckets[category] = genes

    dr_dict = {
        bucket: merged.loc[merged["Systematic ID"].isin(genes), "DR"].dropna().tolist()
        for bucket, genes in buckets.items()
    }

    detail_frames = []
    for bucket, genes in buckets.items():
        sub = final_merged[final_merged["Systematic ID"].isin(genes)].copy()
        sub["Verification result bucket"] = bucket
        detail_frames.append(sub)
    detail = pd.concat(detail_frames, ignore_index=True) if detail_frames else final_merged.iloc[0:0].copy()

    return dr_dict, detail


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


def _boxplot_figure(dr_dict: dict[str, list[float]], title: str) -> plt.Figure:
    """Two-panel figure: boxplot+violin on the left, per-bucket Q1/median/Q3/mean text on the right.

    Byte-faithful to the notebook's boxplot_with_violinplot_and_statistics.
    Empty buckets are dropped so the violin/box call never sees a zero-length
    sample (which matplotlib rejects).
    """
    dr_dict = {k: v for k, v in dr_dict.items() if len(v) > 0}
    fig, axes = plt.subplots(
        1, 2, figsize=(2 * AX_WIDTH, AX_HEIGHT), sharey=True, gridspec_kw={"width_ratios": [3, 1]}
    )
    colors = [CATEGORY_COLOR_MAP.get(bucket, "gray") for bucket in dr_dict]
    boxplot_with_violinplot(list(dr_dict.keys()), list(dr_dict.values()), axes[0], colors)
    axes[0].set_title(f"{title}\nDepletion Rate (DR)")
    axes[0].set_xlim(-0.3, 1.5)

    for row, (bucket, drs) in enumerate(dr_dict.items()):
        q1, median, q3 = np.percentile(drs, [25, 50, 75])
        info = f"Q1={q1:.2f}, Median={median:.2f}, Q3={q3:.2f}, Mean={np.mean(drs):.2f}"
        axes[1].text(0.0, row, info, va="center", ha="left", fontweight="bold")
    axes[1].axis("off")
    fig.tight_layout()
    return fig


def _verification_donut_figure(dr_dict: dict[str, list[float]], title: str) -> plt.Figure:
    """Donut of verification-bucket composition (bucket size = wedge), ordered by _VERIFICATION_BUCKET_ORDER."""
    ordered = {k: dr_dict[k] for k in _VERIFICATION_BUCKET_ORDER if k in dr_dict and len(dr_dict[k]) > 0}
    labels = list(ordered.keys())
    values = [len(v) for v in ordered.values()]
    colors = [DONUT_COLOR_MAP.get(label, "gray") for label in labels]

    fig, ax = plt.subplots(figsize=(AX_WIDTH, AX_HEIGHT))
    if values:
        donut_chart(values=values, labels=labels, colors=colors, center_text="have been\nverified", ax=ax)
    ax.set_title(f"{title}\nverification results")
    fig.tight_layout()
    return fig


def build_boxplot_pdf(
    merged: pd.DataFrame,
    final_merged: pd.DataFrame,
    simplified_verification: pd.DataFrame,
    output_boxplots: Path,
    output_critical_genes_dir: Path,
) -> dict[str, list[str]]:
    """Write the boxplot/violin + donut PDF and per-critical-group review TSVs.

    Returns {group: [outlier Systematic IDs]} so the depletion-curve builder
    plots exactly the same genes each group's boxplot/TSV covers.
    """
    figures: list[plt.Figure] = []

    # Basic per-category boxplot (notebook §4.1): DR grouped by
    # Category_with_essentiality, restricted to the canonical categories.
    basic = merged[merged["cat_canon"].isin(_BASIC_BOXPLOT_CATEGORIES)]
    basic_dict = basic.groupby("Category_with_essentiality")["DR"].apply(list).to_dict()
    figures.append(_boxplot_figure(basic_dict, "Deletion library categories"))

    # Critical-gene groups (notebook §4.2-4.4): boxplot + donut + review TSV each.
    group_genes: dict[str, list[str]] = {}
    for group, spec in _CRITICAL_GROUPS.items():
        dr_dict, detail = prepare_verification_data(
            merged, final_merged, simplified_verification,
            outlier_filter=spec["filter"], sort=spec["sort"],
        )
        group_genes[group] = (
            merged.query(spec["filter"], engine="python")
            .sort_values("DR", ascending=spec["sort"] == "asc")["Systematic ID"]
            .unique()
            .tolist()
        )
        figures.append(_boxplot_figure(dr_dict, group))
        figures.append(_verification_donut_figure(dr_dict, group))
        detail.to_csv(output_critical_genes_dir / f"critical_genes_{group}.tsv", sep="\t", index=False)

    with PdfPages(output_boxplots) as pdf:
        for fig in figures:
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)
    return group_genes


def load_grna_timepoints(grna_path: Path | None) -> pd.DataFrame | None:
    """Load the gRNA per-timepoint LFC table, indexed by Systematic ID (extracted from gRNA_ID).

    gRNA rows are keyed by gene Name in the source file, which is NOT unique;
    the systematic ID embedded in gRNA_ID ("SPAC1002.02_42" -> "SPAC1002.02")
    is unique, so we index on it to align with the DIT-HAP table. Returns None
    when no gRNA file is provided (non-HD datasets render DIT-HAP-only curves).
    """
    if grna_path is None:
        return None
    grna = pd.read_csv(grna_path, index_col=0)
    grna["Systematic ID"] = grna["gRNA_ID"].str.rsplit("_", n=1).str[0]
    return grna.drop_duplicates("Systematic ID").set_index("Systematic ID")


def build_depletion_curve_pdf(
    group_genes: dict[str, list[str]],
    gene_timepoints: pd.DataFrame,
    grna_timepoints: pd.DataFrame | None,
    output_depletion_curves: Path,
) -> None:
    """One 4-column grid per critical group; each panel is a gene's DIT-HAP (+gRNA) depletion curve.

    gene_timepoints must be indexed by Systematic ID and carry A/DR/DL + YES0-4.
    Panels are titled by gene Name (readable), aligned on Systematic ID.
    """
    n_cols = 4
    with PdfPages(output_depletion_curves) as pdf:
        for group, genes in group_genes.items():
            present = [g for g in genes if g in gene_timepoints.index]
            if not present:
                continue
            n_rows = int(np.ceil(len(present) / n_cols))
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * AX_WIDTH, n_rows * AX_HEIGHT))
            axes = np.atleast_1d(axes).flatten()
            for ax, gene in zip(axes, present):
                dit_row = gene_timepoints.loc[gene]
                grna_row = (
                    grna_timepoints.loc[gene]
                    if grna_timepoints is not None and gene in grna_timepoints.index
                    else None
                )
                title = f"{dit_row['Name']} ({gene})" if "Name" in dit_row else gene
                plot_gene_depletion_curve(
                    ax, dit_row, grna_row, title,
                    dit_generations=DIT_HAP_GENERATIONS, grna_generations=GRNA_GENERATIONS,
                )
            for extra_ax in axes[len(present):]:
                fig.delaxes(extra_ax)
            fig.suptitle(f"{group} depletion curves (n={len(present)})")
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)


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
    merged = canonicalize_category(merged)
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

    # Notebook §4-5 additions: boxplot/violin + critical-gene donuts + review
    # TSVs, and (when the timepoint inputs are wired) the depletion curves.
    if config.output_boxplots is not None and config.output_critical_genes_dir is not None:
        verification_full = load_essentiality_verification_full(config.essentiality_verification)
        final_merged = build_final_merged(merged, verification_full)
        group_genes = build_boxplot_pdf(
            merged, final_merged, essentiality_verification,
            config.output_boxplots, config.output_critical_genes_dir,
        )
        if config.output_depletion_curves is not None and config.gene_timepoints is not None:
            gene_tp = load_gene_level(config.gene_timepoints).set_index("Systematic ID")
            grna_tp = load_grna_timepoints(config.grna_timepoints)
            build_depletion_curve_pdf(group_genes, gene_tp, grna_tp, config.output_depletion_curves)

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
    # Notebook §4-5 additions (optional so the simplified two-output invocation still works).
    parser.add_argument("--gene-timepoints", type=Path, default=None, help="Gene-level fitting statistics with YES0-4 (for depletion curves)")
    parser.add_argument("--grna-timepoints", type=Path, default=None, help="gRNA per-timepoint LFC CSV (optional; HD-only overlay)")
    parser.add_argument("--output-boxplots", type=Path, default=None, help="Output boxplot/violin + donut PDF")
    parser.add_argument("--output-depletion-curves", type=Path, default=None, help="Output depletion-curve PDF")
    parser.add_argument("--output-critical-genes-dir", type=Path, default=None, help="Output dir for per-group critical-gene TSVs")
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
            gene_timepoints=args.gene_timepoints,
            grna_timepoints=args.grna_timepoints,
            output_boxplots=args.output_boxplots,
            output_depletion_curves=args.output_depletion_curves,
            output_critical_genes_dir=args.output_critical_genes_dir,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
