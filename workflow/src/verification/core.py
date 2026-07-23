"""
Deletion Library Verification — Core Logic
==========================================

Shared constants, loaders, merge/stats functions, critical-gene analysis, and
figure builders for the verification stage. Ported from
DIT_HAP_pipeline/workflow/notebooks/compare_with_deletion_library.ipynb and
factored out of the original single-script port so the stage can be split into
independent Snakemake rules (prepare -> category summary / boxplots /
depletion curves), each re-runnable on its own.

Design doc: docs/plans/2026-07-22-verification-rules-split-design.md.

NOTE: the curated deletion_library_categories.xlsx schema changed after the
source notebook was written — `Updated_Systematic_ID` no longer exists and
`Systematic ID` now holds the current ID directly; merge_deletion_library()
accepts either column name. The curated `Category` values also drifted (the
notebook's `WT` is now `WT-like`, and several compound multi-phenotype labels
were added). Per project decision the verification stage uses those RAW
Category values verbatim everywhere — display text, boxplot/critical grouping,
and the outlier filters all match the literal curated labels; there is no
folding back to the notebook vocabulary. Colors are the one exception: a raw
label with no direct color entry reuses its phenotype family's representative
color via _category_color_key() (color lookup only, never display text).

Usage
-----
    from workflow.src.verification.core import (
        load_gene_level, load_deletion_library, load_essentiality_verification,
        merge_deletion_library, build_final_merged,
        build_boxplot_pdf, build_depletion_curve_pdf,
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
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
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

# Category display order for the donut / DR-scatter, using the RAW curated
# labels verbatim (no folding). Ordered by phenotype-severity progression
# (most-arrested spores -> healthiest WT-like), with each compound label placed
# next to its leading-phenotype family. Categories absent from the data are
# silently skipped; any raw label NOT listed here is appended after the ordered
# ones (nothing is filtered out).
_CATEGORY_ORDER = [
    "spores",
    "spores, germinated",
    "spores, germinated, divided or microcolonies",
    "spores, miscellaneous",
    "germinated",
    "germinated, divided or microcolonies",
    "microcolonies",
    "microcolonies, small colonies",
    "small colonies (E)",
    "E",
    "very small colonies",
    "small colonies",
    "WT-like",
]



# Legacy -> current metric column names, same quirk as coverage.smk's
# compute_coverage_stats.load_gene_level / clustering's candidates.load_and_annotate:
# some releases' gene-level fitting_results.tsv still ship the pre-rename um/lam
# headers instead of DR/DL.
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}

# COLOR-ONLY fallback (never affects display text). CATEGORY_COLOR_MAP /
# DONUT_COLOR_MAP only know the notebook's original single-phenotype vocabulary,
# but the curated file uses the raw labels "WT-like" and several compound
# multi-phenotype labels that have no direct color entry. For color lookup only,
# each such label reuses its phenotype family's representative color; the label
# itself is always shown verbatim.
_CATEGORY_COLOR_ALIASES = {
    "WT-like": "WT",
    # Leading phenotype determines the family color for the compound labels.
    "spores, germinated": "spores",
    "spores, germinated, divided or microcolonies": "spores",
    "spores, miscellaneous": "spores",
    "germinated, divided or microcolonies": "germinated",
    # "microcolonies, small colonies" -> "small colonies": the finer distinction.
    "microcolonies, small colonies": "small colonies",
}


def _category_color_key(category: str) -> str:
    """Map a raw Category value to a CATEGORY_COLOR_MAP/DONUT_COLOR_MAP key (color lookup only)."""
    return _CATEGORY_COLOR_ALIASES.get(category, category)


# Basic-boxplot category selection using the RAW curated labels verbatim.
# Grouping is by Category_with_essentiality after restricting to these
# categories. Compound multi-phenotype labels are intentionally excluded (they
# are not one of the notebook's canonical single-phenotype buckets).
_BASIC_BOXPLOT_CATEGORIES = ["spores", "germinated", "microcolonies", "very small colonies", "small colonies", "WT-like"]

# The four "critical gene" outlier groups (notebook §4.2-4.4). Each filter runs
# against the RAW `Category` column (literal curated labels — no folding), so
# only genes with those exact single-phenotype labels are selected; the compound
# multi-phenotype labels do not enter these analytical groups. `sort` orders the
# outlier gene list by DR: WT->nonWT / small->E look at highest-DR first (desc),
# E->V lowest (asc).
_CRITICAL_GROUPS = {
    "WT2nonWT": {"filter": "Category == 'WT-like' and DR > 0.35", "sort": "desc"},
    "scE2E": {"filter": "Category == 'small colonies' and DR > 0.75 and DeletionLibrary_essentiality == 'E'", "sort": "desc"},
    "sc2E": {"filter": "Category == 'small colonies' and DR > 0.75 and DeletionLibrary_essentiality != 'E'", "sort": "desc"},
    "E2V": {"filter": "Category in ['spores', 'germinated', 'microcolonies'] and DR < 0.35", "sort": "asc"},
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


def load_deletion_library(deletion_library_path: Path) -> pd.DataFrame:
    """Load the curated deletion library xlsx, keeping just the ID + Category columns.

    Handles both schemas: the old file's `Updated_Systematic_ID` key (also used
    by the unit test fixture) and the current file's `Systematic ID` key.
    """
    deletion_library = pd.read_excel(deletion_library_path)
    id_col = "Updated_Systematic_ID" if "Updated_Systematic_ID" in deletion_library.columns else "Systematic ID"
    return deletion_library[[id_col, "Category"]]


def load_essentiality_verification(essentiality_verification_path: Path) -> pd.DataFrame:
    """Load + rename the curated verification table to the notebook's column names.

    Simplifies compound `verification_phenotype` values down to plain "E" (see
    _VERIFICATION_PHENOTYPE_SIMPLIFY), drops rows missing either call, and
    appends the manual gpd1 row (_MANUAL_VERIFICATION_ROWS).
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
    return pd.concat([verification, _MANUAL_VERIFICATION_ROWS], ignore_index=True)


def load_essentiality_verification_full(essentiality_verification_path: Path) -> pd.DataFrame:
    """Load the curated verification table keeping ALL columns (area day3-6, comments, ...).

    Feeds build_final_merged / the per-critical-group review TSVs, which need the
    colony-area measurements load_essentiality_verification drops. The
    `verification_phenotype` values are NOT simplified here (the raw label is
    what a human reviewer wants to see).
    """
    return pd.read_csv(essentiality_verification_path)


# =============================================================================
# MERGE + STATS (unit-tested)
# =============================================================================
def merge_deletion_library(gene_result: pd.DataFrame, deletion_library: pd.DataFrame) -> pd.DataFrame:
    """Left-merge gene-level results with deletion library categories on Systematic ID.

    Accepts either the old schema (`Updated_Systematic_ID` key) or the current
    schema (`Systematic ID` key) for `deletion_library`.
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
    """Append an "(E)" suffix to 'small colonies' when DeletionLibrary_essentiality == 'E'."""
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
    """Count genes per Category_with_essentiality (see apply_category_with_essentiality)."""
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
    """Compare curated 'Verified essentiality' against 'DeletionLibrary_essentiality' where both known."""
    known = merged_with_verification.dropna(subset=["Verified essentiality", "DeletionLibrary_essentiality"])
    match = int((known["Verified essentiality"] == known["DeletionLibrary_essentiality"]).sum())
    return {
        "verified_total": len(known),
        "match": match,
        "mismatch": len(known) - match,
    }


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
# CRITICAL-GENE ANALYSIS (unit-tested)
# =============================================================================
def build_final_merged(merged: pd.DataFrame, verification_full: pd.DataFrame) -> pd.DataFrame:
    """Right-join gene-level+category data with the FULL verification table (area columns kept).

    Reconstructs the notebook's `final_merged`: one row per curated-verification
    gene, carrying DR/DL/FYPOviability/DeletionLibrary_essentiality/Category plus
    the raw verification phenotype/essentiality and colony-area columns. Genes
    verified essential ('E') but missing a day-3 area are zero-filled for the
    area columns, byte-faithful to the notebook (confirmed dead = zero area).
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

    Selects outliers via `outlier_filter` (run against `merged`, which carries
    the raw `Category` column), crosses them with the simplified verification
    table, buckets
    into {"Not verified": [...], <verified category>: [...]} preserving
    _VERIFICATION_BUCKET_ORDER. Each bucket value is member DR values (boxplot);
    bucket size drives the donut. Second return is the per-gene detail frame
    (from final_merged) tagged with its bucket, for the review TSV.
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


def select_group_outliers(merged: pd.DataFrame, group: str) -> list[str]:
    """Return the DR-sorted, deduped outlier Systematic IDs for a critical group (see _CRITICAL_GROUPS).

    Shared by the boxplot builder (review TSVs) and the depletion-curve builder
    (which genes to plot), so both cover exactly the same gene set.
    """
    spec = _CRITICAL_GROUPS[group]
    return (
        merged.query(spec["filter"], engine="python")
        .sort_values("DR", ascending=spec["sort"] == "asc")["Systematic ID"]
        .unique()
        .tolist()
    )


# =============================================================================
# FIGURES — category summary (donut + DR scatter)
# =============================================================================
def plot_category_donut(category_stats: pd.DataFrame) -> plt.Figure:
    """Donut chart of gene counts per deletion-library phenotype category."""
    raw_categories = set(category_stats["category"])
    ordered = [c for c in raw_categories if c in _CATEGORY_ORDER]
    ordered.sort(key=lambda c: _CATEGORY_ORDER.index(c))
    remaining = [c for c in category_stats["category"] if c not in ordered]
    labels = ordered + remaining
    counts_by_label = category_stats.set_index("category")["count"]
    values = [int(counts_by_label[label]) for label in labels]
    colors = [DONUT_COLOR_MAP.get(_category_color_key(label), "gray") for label in labels]

    unmapped = [(label, int(counts_by_label[label])) for label, color in zip(labels, colors) if color == "gray"]
    if unmapped:
        logger.warning(
            f"{sum(n for _, n in unmapped):,} genes plotted as gray (unmapped category): "
            + ", ".join(f"{label!r} (n={n})" for label, n in unmapped)
        )

    fig, ax = plt.subplots(figsize=(AX_WIDTH, AX_HEIGHT))
    donut_chart(
        values=values, labels=labels, colors=colors,
        center_text=f"Total\n{sum(values):,}\ngenes", ax=ax,
    )
    ax.set_title("Deletion library phenotype categories")
    fig.tight_layout()
    return fig


def plot_dr_scatter_by_category(merged: pd.DataFrame) -> plt.Figure:
    """Scatter of DR per gene, grouped by RAW category (x-jittered), n= on the x-tick labels.

    Every curated category is plotted verbatim — nothing is filtered out. Order
    follows _CATEGORY_ORDER, with any label not listed there appended after the
    ordered ones. Colors use the phenotype-family representative color
    (_category_color_key); a label with no family mapping falls back to gray.
    """
    all_categories = list(merged["Category"].dropna().unique())
    ordered = [c for c in all_categories if c in _CATEGORY_ORDER]
    ordered.sort(key=lambda c: _CATEGORY_ORDER.index(c))
    categories = ordered + [c for c in all_categories if c not in _CATEGORY_ORDER]

    fig, ax = plt.subplots(figsize=(AX_WIDTH, AX_HEIGHT))
    rng = np.random.default_rng(42)
    tick_labels = []
    for i, category in enumerate(categories):
        dr_values = merged.query("Category == @category")["DR"].dropna()
        jitter = rng.uniform(-0.15, 0.15, size=len(dr_values))
        ax.scatter(
            i + jitter, dr_values, alpha=0.5, s=10,
            color=CATEGORY_COLOR_MAP.get(_category_color_key(category), "gray"),
        )
        tick_labels.append(f"{category}\n(n={len(dr_values)})")
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(tick_labels, rotation=30, ha="right")
    ax.set_ylabel("Depletion Rate (DR)")
    ax.set_title("DR by deletion library category")
    fig.tight_layout()
    return fig


def build_category_summary_pdf(category_stats: pd.DataFrame, merged: pd.DataFrame, output_figures: Path) -> None:
    """Write the two-page category-summary PDF: phenotype donut + DR-by-category scatter."""
    fig_donut = plot_category_donut(category_stats)
    fig_scatter = plot_dr_scatter_by_category(merged)
    with PdfPages(output_figures) as pdf:
        pdf.savefig(fig_donut, dpi=300, bbox_inches="tight")
        pdf.savefig(fig_scatter, dpi=300, bbox_inches="tight")
    plt.close(fig_donut)
    plt.close(fig_scatter)


# =============================================================================
# FIGURES — boxplot/violin + critical-group donuts + review TSVs
# =============================================================================
def _boxplot_figure(dr_dict: dict[str, list[float]], title: str) -> plt.Figure:
    """Two-panel figure: boxplot+violin left, per-bucket Q1/median/Q3/mean text right.

    Empty buckets are dropped so the violin/box call never sees a zero-length
    sample (which matplotlib rejects).
    """
    dr_dict = {k: v for k, v in dr_dict.items() if len(v) > 0}
    # Sort by reverse _CATEGORY_ORDER (healthiest -> most arrested)
    _cat_order_index = {cat: i for i, cat in enumerate(_CATEGORY_ORDER)}
    dr_dict = dict(
        sorted(dr_dict.items(), key=lambda item: _cat_order_index.get(item[0], len(_CATEGORY_ORDER)), reverse=True)
    )
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
) -> None:
    """Write the boxplot/violin + critical-group donut PDF and per-group review TSVs."""
    figures: list[plt.Figure] = []

    # Basic per-category boxplot (notebook §4.1): DR grouped by
    # Category_with_essentiality, restricted to the canonical categories.
    basic = merged[merged["Category"].isin(_BASIC_BOXPLOT_CATEGORIES)]
    basic_dict = basic.groupby("Category_with_essentiality")["DR"].apply(list).to_dict()
    figures.append(_boxplot_figure(basic_dict, "Deletion library categories"))

    # Critical-gene groups (notebook §4.2-4.4): boxplot + donut + review TSV each.
    for group, spec in _CRITICAL_GROUPS.items():
        dr_dict, detail = prepare_verification_data(
            merged, final_merged, simplified_verification,
            outlier_filter=spec["filter"], sort=spec["sort"],
        )
        figures.append(_boxplot_figure(dr_dict, group))
        figures.append(_verification_donut_figure(dr_dict, group))
        detail.to_csv(output_critical_genes_dir / f"critical_genes_{group}.tsv", sep="\t", index=False)

    with PdfPages(output_boxplots) as pdf:
        for fig in figures:
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)


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
    merged: pd.DataFrame,
    gene_timepoints: pd.DataFrame,
    grna_timepoints: pd.DataFrame | None,
    output_depletion_curves: Path,
) -> None:
    """One 4-column grid per critical group; each panel a gene's DIT-HAP (+gRNA) depletion curve.

    Genes per group come from select_group_outliers (same set as the boxplots /
    review TSVs). gene_timepoints must be indexed by Systematic ID and carry
    A/DR/DL + YES0-4; panels titled by gene Name, aligned on Systematic ID.
    """
    n_cols = 4
    with PdfPages(output_depletion_curves) as pdf:
        for group in _CRITICAL_GROUPS:
            genes = select_group_outliers(merged, group)
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





