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
BH-correct across all pairs. The figure stage renders one n x n pairwise scatter
matrix (lower triangle, variables clustered into a similarity order) through the
house render_pairwise_matrix_figure and one two-panel correlation heatmap
(Pearson | Spearman) through cns.heatmapplot, both in that same clustered order --
the heatmap is handed the very tree it is ordered by, so it can draw it.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-10-08
Version:  4.1.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import sys
import warnings
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from scipy.stats import false_discovery_control, pearsonr, spearmanr

import matplotlib

matplotlib.use("Agg")  # headless: builders only write figures, never display
import matplotlib.pyplot as plt  # noqa: E402
from loguru import logger  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# Project path setup: sibling src/ modules (figures, figure_render.*) import by
# bare name, which needs workflow/src itself on sys.path.
_SRC_DIR = Path(__file__).resolve().parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.append(str(_SRC_DIR))

from figures import (  # noqa: E402
    apply_house_style,
    house_colors,
    panel_labels,
)

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

# The two correlation coefficients, as (stats column, panel title) pairs: one
# figure holds both, since they are the same matrix read two ways.
HEATMAP_PANELS: tuple[tuple[str, str], ...] = (("r_pearson", "Pearson"), ("rho_spearman", "Spearman"))

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

# Study categories in the order they are keyed in the heatmap legend, and the
# layout of that legend's own cell -- a value colourbar under, category key over,
# both drawn by this module because cns.heatmapplot's own legends hang off the
# last panel and get clipped by their axes (see plot_correlation_heatmap).
STUDY_CATEGORIES = sorted(set(COLUMN_CATEGORIES.values()))

# Correlation is bounded by -1/+1, so the diverging scale is fixed rather than
# fitted to the data: a colour then means the same strength on every panel and
# every run. RdBu_r is the house read of a signed scale (blue negative, red
# positive); the sequential house map is for magnitudes, which this is not.
CORRELATION_MIN = -1
CORRELATION_MAX = 1
HEATMAP_CMAP = "RdBu_r"
# This figure sets its page directly instead of deriving it from a PanelShape,
# because the heatmap is a composite plotter: it lays its body -- matrix, row
# dendrogram, annotation bands, row names -- out against the page and ignores the
# axes box it was handed. A fitted grid_axes panel therefore does not hold a
# house-sized panel, it holds the whole composite, and pushes the page past the
# journal width. Measured with no panel shape in the way, the matrix body comes
# out 97 x 100 -- the house panel, within the 3 px the left annotation band takes
# from the matrix -- at this page, 2 columns, and the rect below. The page is at
# the journal cap: any wider and it is scaled down at typesetting. Change the
# column count or the rect and it has to be re-measured.
HEATMAP_PAGE_PX = (540, 136)

# A white annotation track padded between the Study band and the matrix: the
# plotter packs its annotations straight against the heatmap and has no gap
# setting of its own, so the gap is one more (blank) annotation. The plotter
# labels every track with its column name and offers no way to switch that off,
# hence the blank name and value.
SPACER_COLUMN = " "
SPACER_VALUE = " "
BACKGROUND_COLOR = "#FFFFFF"

# Both legends live in a strip down the right of the page, not inside a panel, so
# neither costs the matrices any width. The strip is reserved with a tight_layout
# rect and the legend axes is placed in that figure-fraction space.
HEATMAP_LEGEND_STRIP_LEFT = 0.84
HEATMAP_LAYOUT_RECT = (0, 0, HEATMAP_LEGEND_STRIP_LEFT, 1)
HEATMAP_LEGEND_BOUNDS = (0.85, 0.02, 0.14, 0.96)

# The category key sits at the top of the strip and the value bar under it, thin
# and upright: at this width a horizontal bar would have to be either short or
# wider than the strip. Legend-axes fraction (x0, y0, w, h).
HEATMAP_CBAR_BOUNDS = (0.22, 0.02, 0.07, 0.40)


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
def correlation_matrix(stats: pd.DataFrame, columns: list[str], metric: str) -> pd.DataFrame:
    """Build the symmetric n x n coefficient matrix of ``columns`` from the long-form stats.

    The diagonal is 1 (a variable against itself) and a pair absent from
    ``stats`` -- one dropped by the overlap/degenerate guards -- keeps its 0
    ("no measured correlation"), which is how the heatmap shows it as a blank
    cell and how the clustering treats it as unrelated.
    """
    matrix = pd.DataFrame(np.eye(len(columns)), index=columns, columns=columns)
    for _, row in stats.iterrows():
        matrix.loc[row["col_x"], row["col_y"]] = row[metric]
        matrix.loc[row["col_y"], row["col_x"]] = row[metric]
    return matrix


def correlation_linkage(stats: pd.DataFrame, columns: list[str], metric: str = "r_pearson"):
    """Average-linkage tree over the 1 - r distances between ``columns``."""
    distance = squareform(np.clip(1 - correlation_matrix(stats, columns, metric).to_numpy(), 0, None), checks=False)
    return linkage(distance, method="average")


def cluster_column_order(
    stats: pd.DataFrame, columns: list[str], metric: str = "r_pearson"
) -> list[str]:
    """Order ``columns`` so similar ones sit together, from the pairwise coefficients.

    Both the matrix figure and the heatmap take this order, so the two figures
    read the same way instead of each clustering on its own (Pearson and Spearman
    rank the columns slightly differently, which would otherwise put the panels
    in different places). The heatmap draws the same tree it is ordered by.
    """
    return [columns[index] for index in leaves_list(correlation_linkage(stats, columns, metric))]


def plot_pairwise_scatter(
    fitness_table: pd.DataFrame,
    columns: list[str],
    output_stem: Path | str,
    *,
    order: list[str],
) -> None:
    """Render the n x n lower-triangle scatter matrix of every pair of ``columns``.

    One figure, not one page per pair: at house panel size the page grows with n
    so every panel stays legible. ``order`` is the clustered column order, which
    also places the companion heatmap, so both figures read consistently.
    """
    from figure_render.scatter import render_pairwise_matrix_figure

    render_pairwise_matrix_figure(
        fitness_table, output_stem, columns=columns, labels=COLUMN_DISPLAY_NAMES, order=order,
    )


def _frame_heatmap_panel(plotter: Any, name: str, letter: str) -> None:
    """Drop every spine of one composite panel, and label it above its own content box.

    The plotter frames the heatmap body with the axes spines and frames one side of
    its annotation bands, which reads as a border that stops halfway. With them off
    the colour cells are the panel. The name and letter cannot go on an axes either:
    the plotter lays its axes out over the same area, so an axes title ends up under
    the body, and the panel's left edge is too close to the page edge for
    ``cns.add_panel_label``'s right-aligned offset. Both are drawn in figure space
    from the panel's own bounding box instead.
    """
    import cnsplots as cns

    figure = plotter.ax.figure
    # Every axes in the figure at this point is the composite's (the legend strip is
    # added afterwards), including the four empty placeholders it keeps for the
    # dendrograms and annotations it is not drawing -- the border is spread over all
    # of them, so the spines are dropped figure-wide rather than by name.
    for member in figure.axes:
        for side in ("left", "right", "top", "bottom"):
            member.spines[side].set_visible(False)

    boxes = [
        member.get_position()
        for member in (plotter.ax_heatmap, plotter.ax_row_dendrogram, plotter.ax_left_annotation)
        if member is not None
    ]
    left, right = min(box.x0 for box in boxes), max(box.x1 for box in boxes)
    top = max(box.y1 for box in boxes)
    above = top + cns.settings.panel_pad_top / (figure.get_size_inches()[1] * 72)

    figure.text(
        left, above, letter,
        ha="left", va="bottom",
        fontsize=cns.settings.title_fontsize,
        fontweight=cns.settings.panel_label_fontweight,
    )
    figure.text(
        (left + right) / 2, above, name,
        ha="center", va="bottom",
        fontsize=cns.settings.title_fontsize,
        fontweight=plt.rcParams["axes.titleweight"],
    )


def _draw_heatmap_legend(cell: Axes, study_colors: dict[str, str], *, fig: Figure) -> None:
    """Draw the shared value colourbar and Study category key inside ``cell`` (the legend strip)."""
    import cnsplots as cns

    cell.set_axis_off()

    colorbar_axes = cell.inset_axes(HEATMAP_CBAR_BOUNDS)
    colorbar = fig.colorbar(
        ScalarMappable(norm=Normalize(vmin=CORRELATION_MIN, vmax=CORRELATION_MAX), cmap=HEATMAP_CMAP),
        cax=colorbar_axes,
        ticks=[CORRELATION_MIN, 0, CORRELATION_MAX],
    )
    colorbar.outline.set_linewidth(cns.settings.axes_linewidth)
    colorbar.set_label("Correlation", labelpad=cns.settings.axes_labelpad)

    handles = [
        Patch(facecolor=color, edgecolor="none", label=category)
        for category, color in study_colors.items()
    ]
    legend = cell.legend(
        handles=handles,
        title="Study",
        loc="upper left",
        bbox_to_anchor=(HEATMAP_CBAR_BOUNDS[0], 1.0),
        frameon=cns.settings.legend_frameon,
        alignment="left",
    )
    legend.get_title().set_ha("left")
    legend.get_title().set_position((0, 0))


def plot_correlation_heatmap(
    stats: pd.DataFrame,
    columns: list[str],
    output_path: Path | str,
    *,
    order: list[str],
) -> None:
    """Render Pearson and Spearman as two panels of one figure, with one shared legend.

    Both coefficients are the same matrix read two ways, so they belong side by
    side and share one value colourbar and one Study category key. Each matrix is
    the house panel size, square, with the row dendrogram to its left and the
    column names dropped (they repeat the row names of the same symmetric matrix).

    The page is chosen rather than derived -- see ``HEATMAP_PAGE_PX`` -- and the
    legends go in a strip down its right edge.

    Both legends are drawn here rather than by cns.heatmapplot: the plotter's own
    legends hang off the last panel's right edge of its axes and get clipped by
    them, whichever page width they are given. Drawing them means owning the
    category colours, which is why they are passed to the heatmap explicitly.
    """
    import anndata as ad
    import cnsplots as cns

    apply_house_style()

    order_labels = [COLUMN_DISPLAY_NAMES.get(c, c) for c in order]
    study_colors = dict(zip(STUDY_CATEGORIES, house_colors(range(len(STUDY_CATEGORIES))), strict=True))
    categories = pd.Categorical(
        [COLUMN_CATEGORIES.get(c, "Other") for c in order],
        categories=STUDY_CATEGORIES,
    )
    # A white annotation track is what separates the Study band from the matrix:
    # the plotter packs its annotations straight against the heatmap and has no
    # gap setting of its own. Row annotation, so the band runs down the side of
    # the matrix where the row labels are, not across the top.
    spacer = pd.Categorical([SPACER_VALUE] * len(order), categories=[SPACER_VALUE])
    annotation_colors = {"Study": study_colors, SPACER_COLUMN: {SPACER_VALUE: BACKGROUND_COLOR}}

    cns.figure(width=HEATMAP_PAGE_PX[0], height=HEATMAP_PAGE_PX[1])
    figure = plt.gcf()
    panel_axes = figure.subplots(1, len(HEATMAP_PANELS), squeeze=False)[0]
    plotters = []

    for ax, (metric, name) in zip(panel_axes, HEATMAP_PANELS, strict=True):
        panel_matrix = correlation_matrix(stats, columns, metric).rename(
            index=COLUMN_DISPLAY_NAMES, columns=COLUMN_DISPLAY_NAMES
        ).loc[order_labels, order_labels]
        data = ad.AnnData(
            panel_matrix.to_numpy(),
            obs=pd.DataFrame({"Study": categories, SPACER_COLUMN: spacer}, index=order_labels),
            var=pd.DataFrame(index=order_labels),
        )
        plotters.append(
            cns.heatmapplot(
                data,
                row_annotation=["Study", SPACER_COLUMN],
                # Clustering is the caller's (both panels are then in one order, and it
                # is the tree the matrix figure is ordered by too), so the tree is
                # handed to the plotter rather than computed per coefficient.
                row_cluster=True,
                row_dendrogram=True,
                row_dendrogram_kws={"linkage": correlation_linkage(stats, columns)},
                col_cluster=False,
                cmap=HEATMAP_CMAP,
                vmin=CORRELATION_MIN,
                vmax=CORRELATION_MAX,
                label="Correlation",
                xlabel="",
                ylabel="",
                xticklabels_rotation=45,
                show_rownames=True,
                show_colnames=False,
                # Explicit colours so the annotation bands match the key drawn below
                # (both panels share the scale, so neither carries its own legend).
                colors=annotation_colors,
                plot_legend=False,
                ax=ax,
            )
        )

    # Lay the composite out inside the part of the page that is not the legend strip.
    figure.tight_layout(rect=HEATMAP_LAYOUT_RECT)
    # Each composite places its own axes during the draw, so the panel labels can
    # only be measured from the boxes they end up in.
    figure.canvas.draw()
    labels = panel_labels(len(HEATMAP_PANELS))
    for plotter, label, (_, name) in zip(plotters, labels, HEATMAP_PANELS, strict=True):
        _frame_heatmap_panel(plotter, name, label)

    _draw_heatmap_legend(figure.add_axes(HEATMAP_LEGEND_BOUNDS), study_colors, fig=figure)
    # The legends are their own axes, which the default tight export bbox drops --
    # saving the full canvas is what keeps them in the figure.
    with cns.settings.context(savefig_bbox="standard"):
        cns.savefig(str(output_path))
    plt.close(plt.gcf())


def plot_comparison_figures(
    fitness_table: pd.DataFrame,
    stats: pd.DataFrame,
    columns: list[str],
    output_dir: Path | str,
) -> None:
    """Render the pairwise matrix figure + the two-coefficient heatmap into ``output_dir``."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    order = cluster_column_order(stats, columns)

    # save_dual() appends .pdf/.review.png itself, so pass the bare stem.
    plot_pairwise_scatter(fitness_table, columns, output_dir / "pairwise_fitness_comparison", order=order)
    plot_correlation_heatmap(stats, columns, output_dir / "correlation_heatmap.pdf", order=order)
