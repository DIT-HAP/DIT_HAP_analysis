"""
Pairwise Fitness Comparison — Core Logic
========================================

Shared constants, loaders, merge/stats functions, and figure builders for the
large-scale-study comparison stage.

Data sources (the "annotated + curated inputs" contract):

- gene_annotation_reference.protein.parquet (1c_annotate) — the DIT-HAP metric
  ``HD_DIT_HAP_DR`` and the gRNA metric ``gRNA_DR`` (sign-flipped to the DIT-HAP
  convention by the annotation stage itself, so no flip happens here).
- pombe_coding_gene_protein_features.tsv (1b_features merge) — the other
  large-scale study fitness/depletion columns (STUDY_FITNESS_COLUMNS).

Pipeline: build a fitness table (features spine + the two merged metrics,
integration-density columns clipped), then correlate every unordered pair of
available fitness columns with Pearson AND Spearman (the density columns stay
heavy-tailed even after clipping, so rank correlation is the robust view) and
BH-correct across all pairs. The figure stage renders the pairwise scatter
matrix through the house render_scatter_grid_figure and a clustered correlation
heatmap through cns.heatmapplot, so both figures read the same stats TSV the
tables do.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-10-08
Version:  3.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import sys
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, pearsonr, spearmanr

import matplotlib

matplotlib.use("Agg")  # headless: builders only write figures, never display
import matplotlib.pyplot as plt  # noqa: E402
from loguru import logger  # noqa: E402

# Project path setup: sibling src/ modules (figures, figure_render.*) import by
# bare name, which needs workflow/src itself on sys.path.
_SRC_DIR = Path(__file__).resolve().parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.append(str(_SRC_DIR))

from figures import apply_house_style  # noqa: E402

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

# The other large-scale study fitness/depletion columns to correlate against,
# byte-faithful to the notebook's fitness_data column list. These live on the
# protein-features table. Column selection at runtime is DEFENSIVE — only those
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

# The two merged metric columns after build_fitness_table, renamed to short
# display names. The DIT-HAP metric comes from the annotation reference's
# gene-level HD_DIT_HAP block; the gRNA metric from its gRNA block (already
# sign-flipped there).
DIT_HAP_FITNESS_COLUMN = "DIT-HAP DR"
GRNA_FITNESS_COLUMN = "gRNA DR"
METRIC_COLUMNS = (DIT_HAP_FITNESS_COLUMN, GRNA_FITNESS_COLUMN)

# Annotation-reference columns the metrics map from.
_SOURCE_METRIC_COLUMNS = {
    "HD_DIT_HAP_DR": DIT_HAP_FITNESS_COLUMN,
    "gRNA_DR": GRNA_FITNESS_COLUMN,
}

# A Pearson/Spearman correlation needs at least this many complete (non-NaN)
# pairs to be meaningful; pairs below this are skipped (logged) rather than
# emitting a degenerate r/p that scipy warns or NaNs on.
MIN_PAIRS_FOR_CORRELATION = 3

# Correlation stats TSV columns in output order.
STATS_COLUMNS = [
    "col_x",
    "col_y",
    "pair",
    "n",
    "r_pearson",
    "p_pearson",
    "p_fdr",
    "rho_spearman",
    "p_spearman",
    "p_spearman_fdr",
]

# Short display names for heatmap tick labels, keyed by full column name.
COLUMN_DISPLAY_NAMES = {
    "Barseq_from_dulab": "Barseq (dulab)",
    "Barseq_from_koch": "Barseq (koch)",
    "Integration density, in-vivo (integrations/kb/million inserts)": "Integration density (in-vivo)",
    "ipkm": "ipkm",
    "uipkm": "uipkm",
    "colony_size_Malecki2016": "Colony size",
    "Max Growth Rate": "Max growth rate",
    "Colony Formation": "Colony formation",
    DIT_HAP_FITNESS_COLUMN: "DIT-HAP DR",
    GRNA_FITNESS_COLUMN: "gRNA DR",
}

# Study category each column belongs to, for the heatmap row annotation bands.
COLUMN_CATEGORIES = {
    "Barseq_from_dulab": "Bar-seq",
    "Barseq_from_koch": "Bar-seq",
    "Integration density, in-vivo (integrations/kb/million inserts)": "Density",
    "ipkm": "Density",
    "uipkm": "Density",
    "colony_size_Malecki2016": "Colony",
    "Max Growth Rate": "Growth",
    "Colony Formation": "Growth",
    DIT_HAP_FITNESS_COLUMN: "This study",
    GRNA_FITNESS_COLUMN: "This study",
}


# Project path setup: src/ modules import their siblings by bare name, which
# needs workflow/src itself on sys.path. Library modules import comparison.core
# as comparison.core only, but the sibling imports (figures, figure_render.*)
# below resolve through this append.
_SRC_DIR = Path(__file__).resolve().parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.append(str(_SRC_DIR))


# =============================================================================
# LOADERS
# =============================================================================
def rename_metrics_for_comparison(annotation_reference: pd.DataFrame) -> pd.DataFrame:
    """Return the annotation reference's two metric columns under short display names.

    Raises if either source column is missing so schema drift surfaces loudly
    instead of silently dropping a metric from the comparison.
    """
    missing = [c for c in _SOURCE_METRIC_COLUMNS if c not in annotation_reference.columns]
    if missing:
        raise KeyError(f"annotation reference is missing metric columns: {missing}")
    return annotation_reference[list(_SOURCE_METRIC_COLUMNS)].rename(
        columns=_SOURCE_METRIC_COLUMNS
    )


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


def compute_correlations(x: pd.Series, y: pd.Series) -> tuple[float, float, float, float]:
    """(Pearson r, Pearson p, Spearman rho, Spearman p) over non-NaN (x, y) pairs.

    Drops any pair with a NaN in either series before calling scipy
    (matches the notebook, which correlates only complete observations).
    """
    paired = pd.DataFrame({"x": x.to_numpy(), "y": y.to_numpy()}).dropna()
    with warnings.catch_warnings():
        # Constant input -> nan stats, guarded by the caller.
        warnings.simplefilter("ignore")
        r, p_value = pearsonr(paired["x"], paired["y"])
        rho, rho_p_value = spearmanr(paired["x"], paired["y"])
    return float(r), float(p_value), float(rho), float(rho_p_value)


# =============================================================================
# CORE LOGIC — merge + assembly
# =============================================================================
def build_fitness_table(
    annotation_reference: pd.DataFrame,
    protein_features: pd.DataFrame,
    clip_upper: float = CLIP_UPPER,
) -> pd.DataFrame:
    """Merge the features spine with the annotation reference's two DR metrics.

    The features table (keyed on gene_systematic_id) is the spine, left-joined
    to the DIT-HAP and gRNA metrics on the same key. The integration-density
    columns are clipped at clip_upper on the way out.
    """
    metrics = rename_metrics_for_comparison(annotation_reference)
    merged = protein_features.merge(
        metrics,
        left_on="gene_systematic_id",
        right_index=True,
        how="left",
    )
    return clip_density_columns(merged, clip_upper=clip_upper)


def select_fitness_columns(fitness_table: pd.DataFrame) -> list[str]:
    """Return the fitness columns present with enough non-NaN data to correlate.

    Defensive against a partial real schema: only STUDY_FITNESS_COLUMNS plus the
    two merged metrics that actually exist in ``fitness_table`` and have at least
    MIN_PAIRS_FOR_CORRELATION non-NaN values are kept. Missing/too-sparse columns
    are logged and skipped so the pairwise loop can't KeyError at runtime.
    """
    candidates = STUDY_FITNESS_COLUMNS + list(METRIC_COLUMNS)
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
    """Long-form Pearson + Spearman stats for every unordered pair of ``columns``.

    Pairs with fewer than MIN_PAIRS_FOR_CORRELATION complete observations are
    skipped (logged) rather than emitting a degenerate r/p. Both families'
    p-values are BH-corrected across the surviving pairs (pearson -> ``p_fdr``,
    spearman -> ``p_spearman_fdr``), so the two FDR columns align row-for-row.
    """
    rows = []
    for col_x, col_y in combinations(columns, 2):
        paired = fitness_table[[col_x, col_y]].dropna()
        n = len(paired)
        if n < MIN_PAIRS_FOR_CORRELATION:
            logger.warning(f"Skipping pair ({col_x} vs {col_y}): only {n} complete pairs")
            continue
        r, p_value, rho, rho_p_value = compute_correlations(paired[col_x], paired[col_y])
        if not all(np.isfinite(v) for v in (r, p_value, rho, rho_p_value)):
            # Constant column(s): the coefficient is undefined, and a NaN p would
            # poison the BH correction downstream.
            logger.warning(f"Skipping pair ({col_x} vs {col_y}): degenerate (constant column)")
            continue
        rows.append({
            "col_x": col_x,
            "col_y": col_y,
            "pair": f"{col_x} vs {col_y}",
            "n": n,
            "r_pearson": r,
            "p_pearson": p_value,
            "rho_spearman": rho,
            "p_spearman": rho_p_value,
        })
    stats = pd.DataFrame(rows, columns=[c for c in STATS_COLUMNS if c not in ("p_fdr", "p_spearman_fdr")])
    if stats.empty:
        return pd.DataFrame(columns=STATS_COLUMNS)

    stats["p_fdr"] = false_discovery_control(stats["p_pearson"].to_numpy(), method="bh")
    stats["p_spearman_fdr"] = false_discovery_control(stats["p_spearman"].to_numpy(), method="bh")
    return stats[STATS_COLUMNS]


# =============================================================================
# PLOTTING
# =============================================================================
def plot_pairwise_scatter(
    fitness_table: pd.DataFrame,
    pairs: list[tuple[str, str]],
    output_stem: Path | str,
) -> None:
    """Render the surviving pairs as a multi-page scatter matrix via the house library.

    ``pairs`` is the list of (col_x, col_y) that SURVIVED the per-pair overlap
    filter in compute_correlation_stats, so the PDF panel set always matches the
    stats TSV row set. Pages are separate PDFs suffixed _p2, _p3, ... (a page
    holds 4 pairs at the house 2-column cap); page 1 keeps the bare stem. Each
    page's df is an inner join of its pairs' columns, so NaN-dropping stays
    local and each panel's n matches its stats TSV row.
    """
    from figure_render.scatter import ScatterPanel, render_scatter_grid_figure

    if not pairs:
        logger.warning("No surviving pairs to plot")
        return

    output_stem = Path(output_stem)
    for page_start in range(0, len(pairs), 4):
        page_number = page_start // 4
        stem = output_stem if page_number == 0 else output_stem.with_name(
            f"{output_stem.stem}_p{page_number + 1}"
        )
        page_pairs = pairs[page_start:page_start + 4]
        page_columns = sorted({c for pair in page_pairs for c in pair})
        page_df = fitness_table[page_columns].dropna()
        panels = [
            ScatterPanel(
                x=col_x,
                y=col_y,
                title="",
                # Axis labels already carry the X-vs-Y reading; a "X vs Y"
                # title collides with the neighbouring panel's on shared pages.
                xlabel=COLUMN_DISPLAY_NAMES.get(col_x, col_x),
                ylabel=COLUMN_DISPLAY_NAMES.get(col_y, col_y),
                show_stats=True,
                density=True,
            )
            for col_x, col_y in page_pairs
        ]
        render_scatter_grid_figure(page_df, stem, panels=panels)


def plot_correlation_heatmap(
    stats: pd.DataFrame,
    columns: list[str],
    metric: str,
    output_path: Path | str,
    *,
    title: str,
) -> None:
    """Clustered square correlation heatmap with category row/col annotation.

    ``metric`` picks the rho/r column of ``stats``; a full symmetric matrix is
    built from the long-form pairs (diagonal = 1). cns.heatmapplot only takes
    category annotations from AnnData, and its lazy imports below
    (anndata/cnsplots) keep the stats-only import path of this module free of
    the plotting stack.
    """
    import anndata as ad
    import cnsplots as cns

    apply_house_style()

    labels = [COLUMN_DISPLAY_NAMES.get(c, c) for c in columns]
    matrix = pd.DataFrame(np.eye(len(columns)), index=labels, columns=labels)
    for _, row in stats.iterrows():
        matrix.loc[
            COLUMN_DISPLAY_NAMES.get(row["col_x"], row["col_x"]),
            COLUMN_DISPLAY_NAMES.get(row["col_y"], row["col_y"]),
        ] = row[metric]
        matrix.loc[
            COLUMN_DISPLAY_NAMES.get(row["col_y"], row["col_y"]),
            COLUMN_DISPLAY_NAMES.get(row["col_x"], row["col_x"]),
        ] = row[metric]

    categories = pd.Categorical(
        [COLUMN_CATEGORIES.get(c, "Other") for c in columns],
        categories=sorted({COLUMN_CATEGORIES[c] for c in columns if c in COLUMN_CATEGORIES}),
    )
    var = pd.DataFrame({"Study": categories}, index=labels)
    data = ad.AnnData(matrix.to_numpy(), obs=pd.DataFrame(index=labels), var=var)

    # ~45 px per row/col of data plus margins for the tick labels; a little
    # savefig padding keeps the annotation legend (drawn at the canvas edge by
    # PyComplexHeatmap) inside the tight bbox instead of clipped.
    side_px = 45 * len(columns) + 230
    cns.figure(width=side_px, height=side_px + 120)
    cns.heatmapplot(
        data,
        col_annotation=["Study"],
        row_cluster=True,
        col_cluster=True,
        cmap="RdBu_r",
        vmin=-1,
        vmax=1,
        label="Correlation",
        xlabel="",
        ylabel="",
        xticklabels_rotation=45,
        show_rownames=True,
        show_colnames=True,
        # Push the annotation legend below the correlation colourbar: both are
        # anchored top-right and would overlap at the default offset.
        legend_vpad=120,
    )
    fig = plt.gcf()
    fig.suptitle(title, y=1.02)
    cns.settings.savefig_pad_inches = 0.4
    try:
        cns.savefig(str(output_path))
    finally:
        cns.settings.savefig_pad_inches = 0.01
    plt.close(fig)


def plot_comparison_figures(
    fitness_table: pd.DataFrame,
    stats: pd.DataFrame,
    columns: list[str],
    output_dir: Path | str,
) -> None:
    """Render the scatter grid + both correlation heatmaps into ``output_dir``."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    surviving_pairs = [
        (str(col_x), str(col_y))
        for col_x, col_y in zip(stats["col_x"], stats["col_y"], strict=True)
    ]
    # save_dual() appends .pdf/.review.png itself, so pass the bare stem.
    plot_pairwise_scatter(fitness_table, surviving_pairs, output_dir / "pairwise_fitness_comparison")
    plot_correlation_heatmap(
        stats, columns, "r_pearson",
        output_dir / "correlation_pearson_heatmap.pdf",
        title="Pearson correlation across large-scale studies",
    )
    plot_correlation_heatmap(
        stats, columns, "rho_spearman",
        output_dir / "correlation_spearman_heatmap.pdf",
        title="Spearman correlation across large-scale studies",
    )
