"""
Verification Figure Renderers
=============================

cnsplots renderers for the deletion-library phenotype verification stage
(``workflow/rules/verification.smk``). Four artifacts, one per analysis:

- ``render_category_summary_figure`` — deletion-library phenotype composition
  (donut) and its DR distribution (strip plot).
- ``render_category_boxplot_figure`` — DR per canonical single-phenotype
  category.
- ``render_critical_group_figure`` — for one critical-gene group: DR per
  verification-result bucket (violin + embedded box) plus that group's
  verification composition (donut).
- ``render_depletion_curves_figure`` — one panel per gene on a single page:
  measured points and the upstream fitted curve for DIT-HAP and, where the
  dataset has one, for the HD gRNA assay overlaid.

Colour contract
---------------
One colour per phenotype, everywhere. The 11 raw curated ``Category`` labels
collapse onto six families (``verification.core.CATEGORY_FAMILY``) because the
house palette holds ten colours and because "spores, germinated" and "spores"
are the same phenotype at different resolution. Families are resolved to Cell
palette entries by index, so a family keeps its colour whichever figure it
appears in; genes with no wet-lab call render in the neutral furniture grey.
cnsplots assigns palette entries positionally, so the panels are painted after
the plot call rather than trusting hue_order.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-17
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import sys
from collections.abc import Sequence
from pathlib import Path

# Add workflow/src to path for imports
_SCRIPT_DIR = Path(__file__).resolve().parent
_SRC_DIR = _SCRIPT_DIR.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
import cnsplots as cns
import matplotlib.pyplot as plt
from loguru import logger
from matplotlib.axes import Axes
from matplotlib.path import Path as MplPath

from figures import (  # noqa: E402
    FURNITURE_COLOR,
    PanelShape,
    apply_house_style,
    fit_panels,
    grid_axes,
    house_colors,
    save_dual,
)
from verification.core import (  # noqa: E402
    BASIC_BOXPLOT_CATEGORIES,
    BUCKET_COLUMN,
    CATEGORY_COLUMN,
    CATEGORY_FAMILIES,
    CATEGORY_WITH_ESSENTIALITY_COLUMN,
    DIT_HAP_FITTED_COLS,
    DIT_HAP_VALUE_COLS,
    DR_COLUMN,
    GENE_NAME_COLUMN,
    GRNA_FITTED_COLS,
    GRNA_RATE_COLUMN,
    GRNA_VALUE_COLS,
    UNVERIFIED_FAMILY,
    category_family,
    order_categories,
)

# =============================================================================
# CONSTANTS
# =============================================================================
# Family -> Cell palette index, in phenotype-severity order (most arrested
# first). Chosen for hue separation across the six, reading warm-to-cool as
# severity falls. Indices are stable positions in HOUSE_PALETTE, not arbitrary
# picks: 0 red, 2 orange, 6 purple, 7 slate, 1 teal, 4 green.
FAMILY_PALETTE_INDICES: dict[str, int] = {
    "essential": 0,
    "spores": 2,
    "germinated": 6,
    "microcolonies": 7,
    "small_colonies": 1,
    "wt_like": 4,
}

# One colour per assay on the depletion-curve panels (Cell palette indices):
# DIT-HAP and the orthologous gRNA. Observed points and the fitted curve share
# their assay's colour and are told apart by marker versus line, so the reader
# does not have to hold four colour meanings in their head.
DIT_HAP_SERIES_PALETTE_INDEX = 0
GRNA_SERIES_PALETTE_INDEX = 1

# Depletion-curve columns per page. 2 columns of WIDE panels is 399 layout px,
# inside the 540 px journal page; the row count follows from the group size,
# since the figure is one page however many genes it holds.
CURVES_PER_ROW = 2

# Fixed y range for every depletion-curve panel, so a gene's depletion can be
# read off one panel and compared against any other. Spans both assays: DIT-HAP
# runs to -9.19 and the gRNA overlay to -14.40 (rmi1, E2V), so a floor of -10 —
# the house default for fitted-curve figures is -10.5 — would clip the deepest
# gRNA curves in all four critical groups. The ceiling is the requested 3; the
# highest value either assay reaches is +1.91.
CURVE_YLIM = (-15.0, 3.0)

# Fraction of cnsplots' embedded-box width to keep. 0.5 halves it, from a third
# of the violin's max width to a sixth. See _narrow_embedded_box.
EMBEDDED_BOX_WIDTH_FACTOR = 0.5

# A scatter point has to stay small enough that the 224-gene WT2nonWT group
# reads as a distribution rather than a smear.
SCATTER_MARKER_SIZE = 8
SCATTER_ALPHA = 0.45

# Category-summary layout. multipanel sizes panels in layout pixels (72/inch)
# from the axes box outwards, so these are the axes boxes plus the space their
# own decorations need. ROW_LABEL_PAD_LEFT_PX is the width of the longest
# label ("spores, germinated, divided or microcolonies  (n=76)", 50 chars at the
# 7 pt house tick size).
FIGURE_WIDTH_PX = 540  # the journal's one-page cap
# Row pitch per panel type. A strip plot only needs to separate its rows; a
# violin has to show a density profile, which the strip pitch is too short for.
STRIP_ROW_HEIGHT_PX = 20
VIOLIN_ROW_HEIGHT_PX = 42
STRIP_PLOT_WIDTH_PX = 120
DR_AXIS_WIDTH_PX = 130
ROW_LABEL_PAD_LEFT_PX = 200
DONUT_PANEL_PX = 150
DONUT_LEGEND_HEIGHT_PX = 30  # room for the legend under the ring
PANEL_LABEL_PAD_PX = 26
FIGURE_DECORATION_PX = 52


def _dr_row_panel(
    panels: cns.multipanel,
    label: str,
    n_rows: int,
    *,
    width: int,
    pad_left: int,
    row_height: int,
) -> Axes:
    """Add a horizontal-DR panel sized to n_rows category rows at the given pitch."""
    return panels.panel(
        label=label,
        width=width,
        height=n_rows * row_height,
        pad_left=pad_left,
        pad_top=PANEL_LABEL_PAD_PX,
        margin_bottom=DONUT_LEGEND_HEIGHT_PX,
    )


# =============================================================================
# COLOR + ORDER HELPERS
# =============================================================================
def family_colors() -> dict[str, str]:
    """Return the family -> hex colour map, resolved from the house palette.

    Read on demand rather than cached at import: ``house_colors`` needs the
    palette that ``apply_house_style`` installs, so a module-level constant
    would capture whatever the palette was before setup ran.
    """
    indices = [FAMILY_PALETTE_INDICES[family] for family in CATEGORY_FAMILIES]
    palette = house_colors(indices)
    colors = dict(zip(CATEGORY_FAMILIES, palette, strict=True))
    # Genes with no wet-lab call are an absence of data, not a phenotype, so
    # they take the neutral furniture grey and stay off the data palette.
    colors[UNVERIFIED_FAMILY] = FURNITURE_COLOR
    return colors


def colors_for_labels(labels: Sequence[str], colors: dict[str, str]) -> list[str]:
    """Return the hex colour for each raw label, logging any label with no family."""
    resolved = []
    for label in labels:
        family = category_family(label)
        if family not in colors:
            logger.warning(f"No colour for label {label!r} (family {family!r}); using grey")
        resolved.append(colors.get(family, FURNITURE_COLOR))
    return resolved


def _paint(artists: Sequence[object], colors: Sequence[str], *, what: str) -> None:
    """Recolour artists positionally, refusing to guess when the counts disagree.

    cnsplots assigns palette entries by position, so a count mismatch means the
    artist->label mapping is not what the caller assumed; painting anyway would
    silently put a phenotype in another phenotype's colour.
    """
    if len(artists) != len(colors):
        logger.warning(f"Expected {len(colors)} {what}, found {len(artists)}; leaving colours as drawn")
        return
    for artist, color in zip(artists, colors, strict=True):
        if hasattr(artist, "set_facecolor"):
            artist.set_facecolor(color)
        if hasattr(artist, "set_color"):
            artist.set_color(color)


def _narrow_embedded_box(ax: Axes) -> None:
    """Shrink each violin's embedded quartile box toward its row centre.

    cnsplots draws the embedded box at 0.2 categorical units against a violin
    max of 0.6, so on these narrow DR distributions the box covers a third of
    the body and hides the density peak it is annotating. The box width is not
    reachable from the public call — ``width=`` is the violin's, and the box's
    lives in a ``boxprops`` dict built inside the function — so the drawn
    artists are scaled here rather than the package's internals duplicated.

    Both the rectangle and the median line are scaled. They are separate
    artists (a patch and a Line2D), so narrowing only the rectangle leaves the
    median overhanging the box it is supposed to sit in.
    """
    def shrink(centre: float, values: np.ndarray) -> np.ndarray:
        return centre + (values - centre) * EMBEDDED_BOX_WIDTH_FACTOR

    for patch in ax.patches:
        path = patch.get_path()
        vertices = path.vertices.copy()
        centre = (vertices[:, 1].min() + vertices[:, 1].max()) / 2
        vertices[:, 1] = shrink(centre, vertices[:, 1])
        patch.set_path(MplPath(vertices, path.codes))

    # A median line is the only box artist with a y extent: whiskers run along
    # the value axis at the row centre, and the violin outline is a collection.
    for line in ax.lines:
        values = np.asarray(line.get_ydata(), dtype=float)
        if values.size and np.ptp(values) > 0:
            centre = (values.min() + values.max()) / 2
            line.set_ydata(shrink(centre, values))


def _paint_legend(ax: Axes, colors: dict[str, str]) -> None:
    """Recolour legend handles to match their label's family colour."""
    legend = ax.get_legend()
    if legend is None:
        return
    for text, handle in zip(legend.get_texts(), legend.legend_handles, strict=True):
        color = colors.get(category_family(text.get_text()))
        if color is not None:
            if hasattr(handle, "set_facecolor"):
                handle.set_facecolor(color)
            if hasattr(handle, "set_color"):
                handle.set_color(color)


def present_in_order(order: Sequence[str], values: Sequence[str] | pd.Series) -> list[str]:
    """Return the labels present in ``values``, ordered by ``order``.

    Anything present but not listed in ``order`` is appended rather than
    dropped: "Not verified" is a bucket of its own but no phenotype, so it has
    no place in the severity order and would otherwise vanish from every panel
    of a critical-group figure.
    """
    present = set(values)
    return [label for label in order if label in present] + sorted(present - set(order))


def parse_generations(value: str) -> list[float]:
    """Parse one ``time_points`` cell (comma-separated generations) into floats."""
    return [float(point) for point in str(value).split(",")]


# =============================================================================
# PANEL RENDERERS
# =============================================================================
def render_gene_depletion_curve_panel(
    ax: Axes,
    dit_row: pd.Series,
    grna_row: pd.Series | None,
    grna_generations: Sequence[float] | None,
    title: str,
) -> None:
    """Draw one gene: both assays' measured points and their upstream fitted curves.

    One colour per assay — DIT-HAP red, gRNA teal — with observed points as
    markers of that colour and the fit as a line of it, so which curve and which
    points belong to which assay is readable without consulting the legend. The
    legend therefore carries only the two fitted curves, each labelled with its
    assay's DR so the two rates compare at a glance; the markers are unlabelled,
    since a second entry per assay would double the legend to say nothing new.

    The curves are the upstream ``*_fitted`` columns, not a re-derivation from
    A/DR/DL. Re-deriving looked equivalent until it met the fits that failed: for
    the 930 genes whose fit is degenerate (DL==0, |fitted| < 0.3) the released
    fit is flat while the model dives to a plateau the data never reaches, so
    roughly a fifth of every group's panels showed a curve contradicting the
    points beside it. The upstream columns are the fit of record.

    ``dit_row`` must carry A/DR/DL, the DIT_HAP_VALUE_COLS, their fitted
    counterparts and a ``time_points`` cell. ``grna_row`` carries the gRNA fit
    parameters and value columns, already sign-flipped into the current
    convention by ``load_grna_timepoints``; it is None for datasets with no gRNA
    assay, which render the DIT-HAP curve alone.
    """
    dit_color = house_colors([DIT_HAP_SERIES_PALETTE_INDEX])[0]
    grna_color = house_colors([GRNA_SERIES_PALETTE_INDEX])[0]

    dit_generations = parse_generations(dit_row["time_points"])
    dit_DR = float(dit_row[DR_COLUMN])
    ax.plot(
        dit_generations,
        dit_row[DIT_HAP_FITTED_COLS].to_numpy(dtype=float),
        color=dit_color,
        label=f"DIT-HAP fit (DR={dit_DR:.3g})",
    )
    # Unlabelled: the markers share their assay's colour with the fitted line, so
    # a legend entry for them would double the legend for no information.
    ax.scatter(
        dit_generations,
        dit_row[DIT_HAP_VALUE_COLS].to_numpy(dtype=float),
        color=dit_color,
        zorder=3,
    )

    if grna_row is not None and grna_generations is not None:
        grna_DR = float(grna_row[GRNA_RATE_COLUMN])
        ax.plot(
            grna_generations,
            grna_row[GRNA_FITTED_COLS].to_numpy(dtype=float),
            color=grna_color,
            label=f"gRNA fit (DR={grna_DR:.3g})",
        )
        ax.scatter(
            grna_generations,
            grna_row[GRNA_VALUE_COLS].to_numpy(dtype=float),
            color=grna_color,
            zorder=3,
        )

    ax.set(title=title, xlabel="Generation", ylabel="LFC", ylim=CURVE_YLIM)
    ax.legend()


def render_dr_boxplot_panel(
    ax: Axes,
    data: pd.DataFrame,
    *,
    category_column: str,
    order: Sequence[str],
    title: str,
    xlabel: str = "Depletion rate (DR)",
) -> None:
    """Violin plot of DR per category with an embedded box, one family colour per violin.

    Horizontal for the same reason as the strip plot: the curated labels are
    long, and on the y axis each gets a row to itself instead of competing for
    the ~20 px of x pitch that seven categories share.

    The violin carries the distribution shape, which a box alone hides at the
    small group sizes here (sc2E buckets hold 3-5 genes).
    """
    cns.violinplot(
        data=data,
        y=category_column,
        x=DR_COLUMN,
        order=list(order),
        # add_box is cnsplots' own embedded quartile box, drawn white inside the
        # coloured body — the package default, and left alone. Overlaying a
        # separate cns.boxplot instead renders the violin with a hole punched
        # through the middle of every body; the cause is the second boxplot call,
        # not the box's colour (a standalone box in any colour reproduces it, and
        # add_box does not).
        add_box=True,
        ax=ax,
    )
    # The violin bodies are the collections; the embedded boxes are patches and
    # keep their own colour.
    _paint(list(ax.collections), colors_for_labels(order, family_colors()), what="violins")
    _narrow_embedded_box(ax)
    # Counts are written here rather than via cns.violinplot's add_count, which
    # formats into the tick labels: this overwrites them, and group size is the
    # one annotation a reader needs to judge a distribution.
    sizes = data[category_column].value_counts()
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{label}  (n={int(sizes.get(label, 0)):,})" for label in order])
    ax.set_ylim(-0.7, len(order) - 0.3)
    # Severity descends down the panel, matching the donut's legend order.
    ax.invert_yaxis()
    ax.set(xlabel=xlabel, ylabel="", title=title)


def render_composition_donut_panel(
    ax: Axes,
    data: pd.DataFrame,
    *,
    category_column: str,
    order: Sequence[str],
    center_text: str,
    title: str,
    legend: str = "bottom",
) -> None:
    """Donut of gene counts per category, wedges in family colours."""
    cns.donutplot(data=data, x=category_column, order=list(order), legend=legend, ax=ax)

    colors = colors_for_labels(order, family_colors())
    _paint(list(ax.patches), colors, what="wedges")

    # cnsplots annotates the hole and the legend title with the column name,
    # which is meaningless for a count donut; replace the former and drop the
    # latter. The panel label may precede it in ax.texts, so match on content.
    for text in ax.texts:
        if text.get_text() == category_column:
            text.set_text(center_text)
    legend_obj = ax.get_legend()
    if legend_obj is not None:
        legend_obj.set_title(None)
    _paint_legend(ax, family_colors())
    ax.set_title(title)


def render_dr_scatter_panel(
    ax: Axes,
    data: pd.DataFrame,
    *,
    category_column: str,
    order: Sequence[str],
    title: str,
    xlabel: str = "Depletion rate (DR)",
) -> None:
    """Strip plot of DR per gene, categories down the y axis, jittered within each row.

    Horizontal because the raw curated labels run to 45 characters ("spores,
    germinated, divided or microcolonies"): on the x axis they either rotate
    into the neighbouring panel or force the panel wider than the page. The
    row labels carry the category, so the panel needs no legend of its own —
    the shared donut legend already maps each category to its colour.
    """
    labels = present_in_order(order, set(data[category_column].dropna()))
    sizes = data.groupby(category_column, observed=True).size()

    rng = np.random.default_rng(42)
    for position, label in enumerate(labels):
        rows = data.loc[data[category_column] == label, DR_COLUMN].dropna()
        jitter = rng.uniform(-0.18, 0.18, size=len(rows))
        ax.scatter(
            rows,
            position + jitter,
            s=SCATTER_MARKER_SIZE,
            alpha=SCATTER_ALPHA,
            color=colors_for_labels([label], family_colors())[0],
        )

    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels([f"{label}  (n={int(sizes.get(label, 0)):,})" for label in labels])
    ax.set_ylim(-0.7, len(labels) - 0.3)
    # Severity descends down the panel, matching the donut legend's order, so a
    # reader can compare the two panels row by row.
    ax.invert_yaxis()
    ax.set(xlabel=xlabel, ylabel="", title=title)


# =============================================================================
# FIGURE ENTRYPOINTS
# =============================================================================
@logger.catch(reraise=True)
def render_category_summary_figure(
    merged: pd.DataFrame,
    output_stem: Path,
    *,
    order: Sequence[str],
) -> None:
    """Write the deletion-library comparison figure: phenotype donut + DR strip plot.

    multipanel rather than a grid: the donut needs a wide, square-ish box for
    its 11-entry legend, while the strip plot is a tall column of 12 rows whose
    row labels alone are ~160 px wide. A uniform grid cannot give one panel a
    wide box and the other a wide left margin.
    """
    apply_house_style()

    n_categories = len(order)
    strip_height_px = n_categories * STRIP_ROW_HEIGHT_PX

    cns.figure(width=FIGURE_WIDTH_PX, height=strip_height_px + FIGURE_DECORATION_PX)
    panels = cns.multipanel(max_width=FIGURE_WIDTH_PX)

    ax_donut = panels.panel(
        label="A",
        width=DONUT_PANEL_PX,
        height=DONUT_PANEL_PX,
        pad_left=PANEL_LABEL_PAD_PX,
        pad_top=PANEL_LABEL_PAD_PX,
        margin_bottom=DONUT_LEGEND_HEIGHT_PX,
    )
    render_composition_donut_panel(
        ax_donut, merged,
        category_column=CATEGORY_COLUMN,
        order=order,
        # The ring only covers genes the curated library classifies; genes
        # absent from the library have no wedge, so the hole must count the
        # ring, not every row of the merged table.
        center_text=f"{int(merged[CATEGORY_COLUMN].notna().sum()):,}\ngenes",
        title="Deletion library phenotype categories",
    )

    ax_strip = _dr_row_panel(
        panels, "B", n_categories,
        width=STRIP_PLOT_WIDTH_PX, pad_left=ROW_LABEL_PAD_LEFT_PX, row_height=STRIP_ROW_HEIGHT_PX,
    )
    render_dr_scatter_panel(
        ax_strip, merged,
        category_column=CATEGORY_COLUMN,
        order=order,
        title="DR by deletion library category",
    )

    save_dual(output_stem)
    logger.success(f"Wrote category summary figure: {output_stem}.pdf")


@logger.catch(reraise=True)
def render_category_boxplot_figure(
    merged: pd.DataFrame,
    output_stem: Path,
) -> None:
    """Write the basic DR violin plot over the canonical single-phenotype categories.

    Restricted to ``BASIC_BOXPLOT_CATEGORIES``: the compound multi-phenotype
    labels are not single buckets, so they do not belong on a categorical axis.
    Splitting "small colonies" by essentiality adds one extra box.
    """
    apply_house_style()

    basic = merged[merged[CATEGORY_COLUMN].isin(BASIC_BOXPLOT_CATEGORIES)]
    boxes = order_categories(basic[CATEGORY_WITH_ESSENTIALITY_COLUMN])

    height_px = len(boxes) * VIOLIN_ROW_HEIGHT_PX
    cns.figure(width=FIGURE_WIDTH_PX, height=height_px + FIGURE_DECORATION_PX)
    panels = cns.multipanel(max_width=FIGURE_WIDTH_PX)
    ax = _dr_row_panel(
        panels, "", len(boxes),
        width=DR_AXIS_WIDTH_PX, pad_left=ROW_LABEL_PAD_LEFT_PX, row_height=VIOLIN_ROW_HEIGHT_PX,
    )
    render_dr_boxplot_panel(
        ax, basic,
        category_column=CATEGORY_WITH_ESSENTIALITY_COLUMN,
        order=boxes,
        title="Deletion library categories",
    )

    save_dual(output_stem)
    logger.success(f"Wrote category box plot figure: {output_stem}.pdf")


@logger.catch(reraise=True)
def render_critical_group_figure(
    dr_by_bucket: dict[str, list[float]],
    output_stem: Path,
    *,
    group: str,
    order: Sequence[str],
) -> None:
    """Write one critical-gene group's figure: DR violin plot + verification donut.

    Buckets that came back empty are dropped before drawing: a zero-length
    sample has no violin and only adds an empty wedge legend entry.
    """
    apply_house_style()

    long = pd.DataFrame(
        [(bucket, value) for bucket, values in dr_by_bucket.items() for value in values],
        columns=[BUCKET_COLUMN, DR_COLUMN],
    )
    present = present_in_order(order, long[BUCKET_COLUMN])
    height_px = max(len(present) * VIOLIN_ROW_HEIGHT_PX, DONUT_PANEL_PX)

    cns.figure(width=FIGURE_WIDTH_PX, height=height_px + FIGURE_DECORATION_PX)
    panels = cns.multipanel(max_width=FIGURE_WIDTH_PX)

    ax_box = _dr_row_panel(
        panels, "A", len(present),
        width=DR_AXIS_WIDTH_PX, pad_left=ROW_LABEL_PAD_LEFT_PX, row_height=VIOLIN_ROW_HEIGHT_PX,
    )
    render_dr_boxplot_panel(
        ax_box, long,
        category_column=BUCKET_COLUMN,
        order=present,
        title=f"{group} — DR by verification result",
    )

    ax_donut = panels.panel(
        label="B",
        width=DONUT_PANEL_PX,
        height=DONUT_PANEL_PX,
        pad_left=PANEL_LABEL_PAD_PX,
        pad_top=PANEL_LABEL_PAD_PX,
        margin_bottom=DONUT_LEGEND_HEIGHT_PX,
    )
    render_composition_donut_panel(
        ax_donut, long,
        category_column=BUCKET_COLUMN,
        order=present,
        center_text=f"{len(long):,}\ngenes\nverified",
        title=f"{group} verification results",
        # Right, not bottom: a critical group has at most a handful of buckets,
        # so the legend fits beside the ring. A bottom legend renders below the
        # axes box, which makes multipanel grow the whole figure and leaves the
        # box plot stranded in a blank band.
        legend="right",
    )

    save_dual(output_stem)
    logger.success(f"Wrote critical group figure for {group}: {output_stem}.pdf")


@logger.catch(reraise=True)
def render_depletion_curves_figure(
    genes: Sequence[str],
    gene_timepoints: pd.DataFrame,
    grna_timepoints: tuple[pd.DataFrame, Sequence[float]] | None,
    output_stem: Path,
    *,
    group: str,
) -> None:
    """Write one critical-gene group's per-gene depletion curves on a single page.

    One page per group, however many genes it holds, so the review PNG covers the
    whole figure rather than its first page. A group can be large — WT2nonWT has
    230 genes — and at two columns that is a very tall page; the alternative
    (extra pages, or a PNG per page) was rejected because the PNG then stops
    being a preview of the artifact.
    """
    apply_house_style()

    present = [gene for gene in genes if gene in gene_timepoints.index]
    if not present:
        logger.warning(f"No genes with per-timepoint data for {group}; nothing to plot")
        return
    if len(present) < len(genes):
        logger.info(f"{len(genes) - len(present)} of {len(genes)} {group} genes have no per-timepoint data")

    grna_rows, grna_generations = grna_timepoints if grna_timepoints is not None else (None, None)

    n_rows = (len(present) + CURVES_PER_ROW - 1) // CURVES_PER_ROW
    logger.info(f"Rendering {len(present)} {group} depletion curves over {n_rows} rows")

    # No panel letters: at this density a label per panel is clutter, and the
    # grid runs past Z anyway.
    axes = grid_axes(n_rows, CURVES_PER_ROW, labels=[], shape=PanelShape.WIDE)

    for ax, gene in zip(axes, present, strict=False):
        dit_row = gene_timepoints.loc[gene]
        grna_row = (
            grna_rows.loc[gene]
            if grna_rows is not None and gene in grna_rows.index
            else None
        )
        name = dit_row[GENE_NAME_COLUMN] if GENE_NAME_COLUMN in dit_row else gene
        # Uncharacterised genes carry their systematic ID as the name, so the
        # parenthetical would just repeat the title.
        title = gene if str(name) == gene else f"{name} ({gene})"
        render_gene_depletion_curve_panel(ax, dit_row, grna_row, grna_generations, title)
    for unused in axes[len(present):]:
        unused.set_visible(False)

    fit_panels()
    save_dual(output_stem)
    logger.success(f"Wrote {n_rows}-row depletion-curve figure: {output_stem}.pdf")
