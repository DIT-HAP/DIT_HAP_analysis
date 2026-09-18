"""
Deletion Library Verification — Core Logic
==========================================

Shared constants, loaders, merge/stats functions and critical-gene selection for
the verification stage. Ported from
DIT_HAP_pipeline/workflow/notebooks/compare_with_deletion_library.ipynb and
factored so the stage can be split into independent Snakemake rules
(prepare -> category summary / boxplots / depletion curves), each re-runnable on
its own.

Data logic only — everything that draws lives in
``workflow/src/figure_render/verification.py`` (ADR-0001: computation and
plotting are separate layers).

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
label resolves to its phenotype family (CATEGORY_FAMILY) for *color* only,
never for display text or grouping.

Usage
-----
    from verification.core import (
        load_deletion_library, load_verification,
        merge_deletion_library, build_final_merged, select_group_outliers,
    )
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from collections.abc import Iterable
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
from release_schema import GRNA_METRIC_SIGN, read_gene_level

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

# Phenotype families, most-arrested first. Every raw curated label and every
# verification-result bucket resolves to exactly one of these (or to
# "unverified" for genes with no wet-lab call). Display text is always the raw
# label; the family exists so one phenotype keeps one colour across every panel
# of the stage (design decision 2026-09-17 — Cell, the house palette, holds 10
# colours and the raw vocabulary has 11 labels).
CATEGORY_FAMILIES = ("essential", "spores", "germinated", "microcolonies", "small_colonies", "wt_like")

# Genes with no wet-lab call are not a phenotype — the renderer draws them in the
# neutral furniture grey rather than a palette colour.
UNVERIFIED_FAMILY = "unverified"

# Raw curated label -> family. Compound multi-phenotype labels follow their
# leading (most severe) phenotype. "small colonies (E)" is the
# Category_with_essentiality spelling of an essential call, so it joins
# "essential" rather than "small_colonies".
CATEGORY_FAMILY = {
    "E": "essential",
    "E (tiny colonies)": "essential",
    "small colonies (E)": "essential",
    "spores": "spores",
    "spores, germinated": "spores",
    "spores, germinated, divided or microcolonies": "spores",
    "spores, miscellaneous": "spores",
    "germinated": "germinated",
    "germinated, divided or microcolonies": "germinated",
    "microcolonies": "microcolonies",
    "microcolonies, small colonies": "microcolonies",
    "small colonies": "small_colonies",
    "very small colonies": "small_colonies",
    "WT": "wt_like",
    "WT-like": "wt_like",
    "Not verified": UNVERIFIED_FAMILY,
}

# Basic-boxplot category selection using the RAW curated labels verbatim.
# Grouping is by Category_with_essentiality after restricting to these
# categories. Compound multi-phenotype labels are intentionally excluded (they
# are not one of the notebook's canonical single-phenotype buckets).
BASIC_BOXPLOT_CATEGORIES = ["spores", "germinated", "microcolonies", "very small colonies", "small colonies", "WT-like"]

# The four "critical gene" outlier groups (notebook §4.2-4.4). Each filter runs
# against the RAW `Category` column (literal curated labels — no folding), so
# only genes with those exact single-phenotype labels are selected; the compound
# multi-phenotype labels do not enter these analytical groups. `sort` orders the
# outlier gene list by DR: WT->nonWT / small->E look at the most depleted first,
# E->V the least (sign flipped 2026-09-17 — negative DR is the depleted end, so
# both the thresholds and every `sort` direction are mirrored from the original).
CRITICAL_GROUPS = {
    "WT2nonWT": {"filter": "Category == 'WT-like' and DR < -0.35", "sort": "asc"},
    "scE2E": {"filter": "Category == 'small colonies' and DR < -0.75 and DeletionLibrary_essentiality == 'E'", "sort": "asc"},
    "sc2E": {"filter": "Category == 'small colonies' and DR < -0.75 and DeletionLibrary_essentiality != 'E'", "sort": "asc"},
    "E2V": {"filter": "Category in ['spores', 'germinated', 'microcolonies'] and DR > -0.35", "sort": "desc"},
}

# Verification-result bucket order for the critical-group boxplots/donuts,
# byte-faithful to the notebook's prepare_verification_data loop + label_orders.
_VERIFICATION_BUCKET_ORDER = [
    "spores", "germinated", "microcolonies", "E", "E (tiny colonies)",
    "very small colonies", "small colonies", "WT",
]

# Bucket shown for critical-group genes that carry no wet-lab call.
NOT_VERIFIED_BUCKET = "Not verified"

# gpd1: the source notebook manually appends this gene to the simplified
# verification table (cell 3) as an E/E verified call — it wasn't in the curated
# verification file then and still isn't now, but it IS a real release gene, so
# the append stays meaningful. Kept byte-faithful.
_MANUAL_VERIFICATION_ROWS = pd.DataFrame(
    {"Systematic ID": ["SPBC215.05"], "Verification result": ["E"], "verification_essentiality": ["E"]}
)

# DIT-HAP per-timepoint LFC columns in the release gene-level statistics table
# (5 time points; the assay always has this shape), and the same points as
# fitted by the upstream curve fit.
DIT_HAP_VALUE_COLS = ["YES0", "YES1", "YES2", "YES3", "YES4"]
DIT_HAP_FITTED_COLS = [f"{column}_fitted" for column in DIT_HAP_VALUE_COLS]

# Column names of the prepared tables, named once so a rename is one edit.
CATEGORY_COLUMN = "Category"
CATEGORY_WITH_ESSENTIALITY_COLUMN = "Category_with_essentiality"
DR_COLUMN = "DR"
GENE_NAME_COLUMN = "Name"
BUCKET_COLUMN = "Verification result bucket"

# gRNA per-timepoint LFC columns in the curated HD gRNA fitted-parameters table
# (6 time points — the gRNA assay samples one more than DIT-HAP).
GRNA_VALUE_COLS = ["M_G0Tet", "M_YES1_Tet", "M_YES2_Tet", "M_YES3_Tet", "M_YES4_Tet", "M_YES5_Tet"]

# The same six time points as fitted by the upstream curve fit.
GRNA_FITTED_COLS = [f"{column}_fitted" for column in GRNA_VALUE_COLS]

# Fit parameters in the curated gRNA table, mirroring the release's A/DR/DL.
GRNA_AMPLITUDE_COLUMN = "A"
GRNA_RATE_COLUMN = "um"
GRNA_LAG_COLUMN = "lam"

# The curated gRNA table is frozen at the pre-2026-09-17 sign convention —
# positive = depleted — while upstream flipped DIT-HAP so negative is now the
# depleted end. Flip the gRNA sign-carrying columns on the way in so both assays
# point the same way on one panel; without it the gRNA curve rises where
# DIT-HAP falls.
#
# The amplitude is flipped with the values, not just the values: in this fitting
# convention A *is* the plateau the curve settles at, so it carries the curve's
# direction. DIT-HAP's A is negative for depleted genes; gRNA's is positive, and
# a curve drawn from the raw A would run away from the DIT-HAP one even with the
# LFC values negated.
#
# ``lam`` is a lag in generations, not a signed magnitude, and is left alone.
GRNA_SIGN_FLIP_COLS = [GRNA_AMPLITUDE_COLUMN, GRNA_RATE_COLUMN, *GRNA_VALUE_COLS, *GRNA_FITTED_COLS]


def category_family(label: str | None) -> str:
    """Resolve a raw curated label or verification bucket to its phenotype family.

    Unknown labels fall back to ``wt_like`` only if literally blank; anything else
    returns itself so an unmapped label shows up as a missing colour rather than
    silently sharing a family's meaning. Callers treat an unknown family as
    "draw in furniture grey" and log it.
    """
    if label is None or (isinstance(label, float) and pd.isna(label)):
        return UNVERIFIED_FAMILY
    return CATEGORY_FAMILY.get(str(label), str(label))


def order_labels_present(known_order: list[str], values: Iterable[str]) -> list[str]:
    """Return the labels present in `values`, ordered by `known_order`, unlisted ones appended.

    Anything present but unlisted is appended rather than dropped: "Not verified" is a
    bucket of its own but no phenotype, so it has no place in the severity order and
    would otherwise vanish from every panel of a critical-group figure.

    `values` may be a Series, array or set — the renderer passes whichever it has.
    """
    present = {v for v in pd.Series(list(values)).dropna()}
    return [label for label in known_order if label in present] + sorted(present - set(known_order))


def order_categories(labels: pd.Series) -> list[str]:
    """Order raw curated Category labels by phenotype severity (see _CATEGORY_ORDER)."""
    return order_labels_present(_CATEGORY_ORDER, labels)


def order_verification_buckets(labels: pd.Series) -> list[str]:
    """Order verification-result buckets from most-arrested to healthiest (see _VERIFICATION_BUCKET_ORDER)."""
    return order_labels_present(_VERIFICATION_BUCKET_ORDER, labels)


# =============================================================================
# LOADERS
# =============================================================================
def load_deletion_library(deletion_library_path: Path) -> pd.DataFrame:
    """Load the curated deletion library xlsx, keeping just the ID + Category columns.

    The pre-2026-07 schema spelled the id column `Updated_Systematic_ID`; it is renamed
    here so every consumer sees one column name (this is the only place the two schemas
    are told apart).
    """
    deletion_library = pd.read_excel(deletion_library_path)
    if "Updated_Systematic_ID" in deletion_library.columns:
        deletion_library = deletion_library.rename(columns={"Updated_Systematic_ID": "Systematic ID"})
    return deletion_library[["Systematic ID", "Category"]]


def load_verification(essentiality_verification_path: Path) -> pd.DataFrame:
    """Load the curated verification table: every column, plus the simplified phenotype.

    One table for both uses the stage makes of it. The per-critical-group review
    TSVs want the raw ``verification_phenotype`` and the colony-area columns; the
    outlier bucketing wants the compound phenotype labels collapsed to plain "E"
    (see _VERIFICATION_PHENOTYPE_SIMPLIFY). Loading it twice to produce two
    projections was the stage's only reason to read the same CSV twice.

    ``systematic_id`` is renamed to the ``Systematic ID`` key the rest of the
    stage joins on; the curated names are otherwise kept as they are, so the
    essentiality column stays ``verification_essentiality``. Rows missing either
    call are dropped, and the notebook's hand-added gpd1 row is appended.
    """
    verification = pd.read_csv(essentiality_verification_path).rename(
        columns={"systematic_id": "Systematic ID"}
    )
    verification["Verification result"] = verification["verification_phenotype"].replace(
        _VERIFICATION_PHENOTYPE_SIMPLIFY
    )
    verification = verification.dropna(subset=["Verification result", "verification_essentiality"])
    return pd.concat([verification, _MANUAL_VERIFICATION_ROWS], ignore_index=True)


def load_grna_timepoints(grna_path: Path | None) -> tuple[pd.DataFrame, list[float]] | None:
    """Load the curated HD gRNA fitted-parameters table, indexed by Systematic ID.

    Source: ``resources/curated/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv``
    — one row per gene (the order-1 gRNA already selected upstream), so no
    de-duplication is needed; the previous source was a per-gRNA export that had
    to be keyed off the systematic ID embedded in ``gRNA_ID``.

    The generation grid is read per row from the ``time_points`` CSV string
    rather than hardcoded: the gRNA assay's six sample times differ from
    DIT-HAP's five, and every row must agree on them for one overlay to be
    meaningful.

    The value, amplitude and rate columns are sign-flipped by
    ``GRNA_METRIC_SIGN`` (see GRNA_SIGN_FLIP_COLS): the curated file is frozen
    at the pre-2026-09-17 convention where positive means depleted, while
    DIT-HAP's columns now run negative for depletion. Without the flip the two
    curves on one panel would rise and fall in opposite directions — see
    comparison/core.py, where the same flip is applied.

    Returns ``(frame, generations)``, or None when no path is given (non-HD
    datasets render DIT-HAP-only curves).
    """
    if grna_path is None:
        return None

    grna = pd.read_csv(grna_path, sep="\t")
    missing = [c for c in ("Systematic ID", "time_points", *GRNA_SIGN_FLIP_COLS) if c not in grna.columns]
    if missing:
        raise KeyError(f"gRNA table {grna_path} is missing columns: {missing}")

    grids = grna["time_points"].astype(str).str.split(",").map(lambda pts: [float(p) for p in pts])
    if grids.map(tuple).nunique() != 1:
        raise ValueError(f"gRNA table {grna_path} carries more than one time-point grid")
    generations = grids.iloc[0]

    if len(generations) != len(GRNA_VALUE_COLS):
        raise ValueError(
            f"gRNA table {grna_path} has {len(generations)} time points but "
            f"{len(GRNA_VALUE_COLS)} value columns ({GRNA_VALUE_COLS})"
        )

    grna[GRNA_SIGN_FLIP_COLS] = grna[GRNA_SIGN_FLIP_COLS] * GRNA_METRIC_SIGN
    indexed = grna.drop_duplicates("Systematic ID").set_index("Systematic ID")
    logger.info(f"Loaded gRNA time points for {len(indexed):,} genes at generations {generations}")
    return indexed, generations


# =============================================================================
# MERGE + STATS (unit-tested)
# =============================================================================
def merge_deletion_library(gene_result: pd.DataFrame, deletion_library: pd.DataFrame) -> pd.DataFrame:
    """Left-merge gene-level results with deletion library categories on Systematic ID.

    `deletion_library` is expected to come from load_deletion_library, which normalizes
    the legacy `Updated_Systematic_ID` spelling.
    """
    return gene_result.merge(deletion_library, on="Systematic ID", how="left")


def apply_category_with_essentiality(row: pd.Series) -> str:
    """Append an "(E)" suffix to 'small colonies' when DeletionLibrary_essentiality == 'E'."""
    if row["Category"] == "small colonies" and row["DeletionLibrary_essentiality"] == "E":
        return f"{row['Category']} (E)"
    return row["Category"]


def count_by(merged: pd.DataFrame, column: str) -> pd.DataFrame:
    """Count genes per value of ``column`` as a category/count frame."""
    return (
        merged.groupby(column, dropna=False)
        .size()
        .reset_index(name="count")
        .rename(columns={column: "category"})
    )


def compute_verification_match_stats(merged: pd.DataFrame, verification: pd.DataFrame) -> dict[str, int]:
    """Compare the curated verification call against the library's, where both are known.

    Only genes carrying both calls can agree or disagree, so this is an inner
    join: a gene the library never classified, or the wet-lab work never reached,
    is not evidence either way.
    """
    both = (
        merged[["Systematic ID", "DeletionLibrary_essentiality"]]
        .merge(verification[["Systematic ID", "verification_essentiality"]], on="Systematic ID", how="inner")
        .dropna()
    )
    match = int((both["verification_essentiality"] == both["DeletionLibrary_essentiality"]).sum())
    return {"verified_total": len(both), "match": match, "mismatch": len(both) - match}


def build_stats_table(merged: pd.DataFrame, verification: pd.DataFrame) -> pd.DataFrame:
    """Flatten category counts + verification match/mismatch counts into one long-form table."""
    rows = [
        {"metric": f"{label}_count", "category": row["category"], "count": row["count"]}
        for label, column in (
            ("category", CATEGORY_COLUMN),
            ("category_with_essentiality", CATEGORY_WITH_ESSENTIALITY_COLUMN),
        )
        for _, row in count_by(merged, column).iterrows()
    ]
    rows += [
        {"metric": "verification", "category": key, "count": count}
        for key, count in compute_verification_match_stats(merged, verification).items()
    ]
    return pd.DataFrame(rows)


# =============================================================================
# CRITICAL-GENE ANALYSIS (unit-tested)
# =============================================================================
def build_final_merged(merged: pd.DataFrame, verification: pd.DataFrame) -> pd.DataFrame:
    """Right-join gene-level+category data with the verification table (area columns kept).

    Reconstructs the notebook's `final_merged`: one row per curated-verification
    gene, carrying DR/DL/FYPOviability/DeletionLibrary_essentiality/Category plus
    the raw verification phenotype/essentiality and colony-area columns. Genes
    verified essential ('E') but missing a day-3 area are zero-filled for the
    area columns, byte-faithful to the notebook (confirmed dead = zero area).
    """
    area_cols = [c for c in verification.columns if "area" in c]
    final = merged.merge(verification, on="Systematic ID", how="right")
    e_missing_day3 = final.query(
        "verification_essentiality == 'E' and median_area_day3.isna()",
        engine="python",
    ).index
    final.loc[e_missing_day3, area_cols] = 0
    return final


def prepare_verification_data(
    merged: pd.DataFrame,
    final_merged: pd.DataFrame,
    verification: pd.DataFrame,
    outlier_filter: str,
    sort: str = "desc",
) -> tuple[dict[str, list[float]], pd.DataFrame]:
    """Bucket a group's outliers by verification result; return {bucket: [DR...]} + gene detail.

    Selects outliers via `outlier_filter` (run against `merged`, which carries
    the raw `Category` column), crosses them with the verification table, buckets
    into {"Not verified": [...], <verified category>: [...]} preserving
    _VERIFICATION_BUCKET_ORDER. Each bucket value is member DR values (boxplot);
    bucket size drives the donut. Second return is the per-gene detail frame
    (from final_merged) tagged with its bucket, for the review TSV.
    """
    outliers = _sorted_unique_ids(merged, outlier_filter, sort == "asc")
    verified = verification[verification["Systematic ID"].isin(outliers)]
    verified_genes = set(verified["Systematic ID"])
    missing = [g for g in outliers if g not in verified_genes]

    buckets: dict[str, list[str]] = {}
    if missing:
        buckets[NOT_VERIFIED_BUCKET] = missing
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
        # Verified buckets come from the verification join. The unverified bucket
        # has no wet-lab record to join to, so its rows come from the gene-level
        # table and carry NaN for every verification column — otherwise the
        # figure advertises "Not verified (n=7)" and the review TSV has no such
        # rows to look at.
        source = merged if bucket == NOT_VERIFIED_BUCKET else final_merged
        sub = source[source["Systematic ID"].isin(genes)].copy()
        sub[BUCKET_COLUMN] = bucket
        detail_frames.append(sub)
    detail = pd.concat(detail_frames, ignore_index=True) if detail_frames else final_merged.iloc[0:0].copy()

    return dr_dict, detail


def _sorted_unique_ids(merged: pd.DataFrame, filter_expr: str, ascending: bool) -> list[str]:
    """The DR-sorted, deduped Systematic IDs matching `filter_expr`.

    One implementation because the boxplot builder (review TSVs) and the depletion-curve
    builder (which genes to plot) must cover exactly the same gene set, and they used to
    hold a copy each — one edit to either filter silently desynchronised them.
    """
    return (
        merged.query(filter_expr, engine="python")
        .sort_values("DR", ascending=ascending)["Systematic ID"]
        .unique()
        .tolist()
    )


def select_group_outliers(merged: pd.DataFrame, group: str) -> list[str]:
    """Return the DR-sorted, deduped outlier Systematic IDs for a critical group (see CRITICAL_GROUPS)."""
    spec = CRITICAL_GROUPS[group]
    return _sorted_unique_ids(merged, spec["filter"], spec["sort"] == "asc")


def critical_group_boxplot_data(
    merged: pd.DataFrame,
    final_merged: pd.DataFrame,
    verification: pd.DataFrame,
    group: str,
) -> tuple[dict[str, list[float]], pd.DataFrame]:
    """``prepare_verification_data`` for one named critical group (see CRITICAL_GROUPS)."""
    spec = CRITICAL_GROUPS[group]
    return prepare_verification_data(
        merged, final_merged, verification,
        outlier_filter=spec["filter"], sort=spec["sort"],
    )
