"""
Shared Rendering Library
========================

House style, panel layout and dual-artifact saving for cnsplots figures.
This is a library module (no main(), no CLI) per python-script-conventions.

Every uniform grid in the project goes through ``grid_axes``, which is the
cnsplots-documented matplotlib grid pattern: ``cns.figure()`` sizes the page,
``fig.subplots()`` allots equal GridSpec cells, ``set_box_aspect`` keeps them
square and ``tight_layout`` fits the decorations. Heterogeneous figures use
``cns.multipanel`` directly.

Sizing runs panel-first: titles, tick labels and axis labels are fixed at the
house font sizes, so legibility depends on the *axes* box, not the page. A
caller picks a ``PanelShape`` and the figure size follows from the grid, rather
than dividing a fixed page width by the column count and letting panels shrink
out of the legible range as the grid grows.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-01
Version:  3.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
from loguru import logger
from matplotlib.axes import Axes


# =============================================================================
# CONSTANTS
# =============================================================================
# Page width in cnsplots layout pixels (72 per inch). Taken from cnsplots' own
# multipanel default. No longer sizes grids -- it is the journal's one-page cap
# that derived widths are checked against.
JOURNAL_WIDTH_PX = int(cns.settings.multipanel_max_width)
JOURNAL_HEIGHT_PX = 425  # ~150 mm, 4:3 aspect


class PanelShape(StrEnum):
    """The three permitted axes-box footprints, in layout pixels."""

    SQUARE = "square"  # 100 x 100
    TALL = "tall"  # 100 wide x 150 high
    WIDE = "wide"  # 150 wide x 100 high


# Axes-box (width, height) per shape. 100-150 px is the legible band at the
# house font sizes: below ~80 px a title or tick label is wider than the box it
# annotates, above ~150 px the type looks undersized against the data area.
PANEL_SIZES_PX: dict[PanelShape, tuple[int, int]] = {
    PanelShape.SQUARE: (100, 100),
    PanelShape.TALL: (100, 150),
    PanelShape.WIDE: (150, 100),
}

# Layout pixels each cell needs on top of its axes box for decorations: y ticks +
# ylabel + panel label across, title + x ticks + xlabel down. Constant across
# grid sizes, since the house font sizes are fixed.
PANEL_DECORATION_WIDTH_PX = 40
PANEL_DECORATION_HEIGHT_PX = 31

# tight_layout's outer pad, which is charged once per figure rather than per cell.
# Solved from measured axes boxes at 1 and 2 columns: the per-cell terms above
# are the slope, this is the intercept. Omitting it made a 1x1 figure's axes
# 89 px instead of 100, since one cell absorbs the whole pad.
FIGURE_PAD_PX = 19

# fit_panels' measure-and-correct loop. Two passes take the worst measured case
# (4 columns, 15 px short) to under a pixel; the tolerance is well below the
# smallest visible difference, and stopping there avoids a pass that would only
# chase rounding.
_FIT_PASSES = 3
_FIT_TOLERANCE_PX = 0.5

# Where grid_axes parks the grid description for fit_panels. On the figure rather
# than in a module global so two figures built in one process cannot cross.
_GRID_SPEC_ATTR = "_house_grid_spec"


@dataclass(kw_only=True, slots=True, frozen=True)
class _GridSpecification:
    """What fit_panels needs to know about the grid grid_axes built."""

    n_rows: int
    n_cols: int
    shape: PanelShape
    axes: tuple[Axes, ...]

# Panel-label offsets from the axes corner, in layout pixels. cnsplots defaults
# both to 0, which puts the letter on top of the y tick labels; measured over a
# 2x6 log grid, this is the smallest pair with no overlap against either the
# tick labels or the axes title. Increased pad_top from 6 to 10 to prevent
# overlap with longer panel titles.
PANEL_LABEL_PAD_LEFT_PX = 25
PANEL_LABEL_PAD_TOP_PX = 25

# Uppercase A-Z, then A1, A2, ... A bare chr(65 + i) emits '[', '\\', ']' past Z.
_ALPHABET_SIZE = 26

# One qualitative palette for every figure in the project. A reader assumes
# consistent colour means the same thing, so the palette name lives here alone
# and no module may name a colour of its own.
HOUSE_PALETTE = "Cell"

# Semantic colour roles, resolved from HOUSE_PALETTE on demand rather than
# written as literals, so re-pointing HOUSE_PALETTE re-colours every figure.
# Indices are stable positions in the palette, not arbitrary picks:
#   0, 1  the first two entries, for a two-series contrast
#   2,3,4 the greyscale-safe trio (min Rec. 601 luma gap 0.145, against 0.026
#         for 0/1/2) — figures get printed and photocopied
_OBSERVED_INDEX = 0
_FITTED_INDEX = 1
_SERIES_INDICES = (2, 3, 4)

# Reference lines, diagonals, cutoff markers, bar edges and annotations are
# furniture, not data. Keeping them off-palette leaves the palette to mean data.
FURNITURE_COLOR = cns.GRAY

# Density is a magnitude, so it needs a sequential map. Installed as
# cns.settings.palette_seq by apply_house_style(), which lands it in
# rcParams['image.cmap'], so unqualified scatter/imshow calls inherit it.
# Overrides cnsplots' 'gnuplot', which is non-monotonic in perceptual lightness:
# equal density steps must read as equal colour steps.
DENSITY_CMAP = "viridis"

# cnsplots' style sets xtick.minor.size == xtick.major.size, so reading
# minor.size back gives ticks indistinguishable from major ones. Derive the
# minor length from major instead, at matplotlib's own default 2.0/3.5 ratio.
_MINOR_TICK_LENGTH_RATIO = 0.57


# =============================================================================
# CORE LOGIC
# =============================================================================
def _style_minor_ticks(ax: Axes) -> None:
    """Shorten minor ticks relative to major ones so the two read as distinct."""
    ax.tick_params(
        which="minor",
        length=plt.rcParams["xtick.major.size"] * _MINOR_TICK_LENGTH_RATIO,
        width=plt.rcParams["xtick.minor.width"],
    )


def apply_log_scale(ax: Axes, *, x: bool = True, y: bool = True) -> None:
    """Put the requested axes on log scale with automatic minor ticks.

    The single place log axes are configured, so every renderer gets the same
    treatment rather than each calling bare ``set_xscale('log')``.

    Minor ticks are styled to be shorter than major ticks. If both x and y are
    scaled and matplotlib's automatic minor tick placement is inconsistent
    (one axis gets minors, the other doesn't), both are disabled for uniformity.
    """
    axes_to_scale = []
    if x:
        ax.set_xscale("log")
        axes_to_scale.append(ax.xaxis)
    if y:
        ax.set_yscale("log")
        axes_to_scale.append(ax.yaxis)

    if not axes_to_scale:
        return

    _style_minor_ticks(ax)

    # If both axes are log-scaled, check for consistency
    if len(axes_to_scale) == 2:
        # Force a draw to populate the minor locators
        ax.figure.canvas.draw()

        x_minors = len(ax.xaxis.get_minor_locator()())
        y_minors = len(ax.yaxis.get_minor_locator()())

        # If inconsistent (one has minors, the other doesn't), disable both
        if (x_minors > 0) != (y_minors > 0):
            ax.xaxis.set_minor_locator(plt.NullLocator())
            ax.yaxis.set_minor_locator(plt.NullLocator())


def apply_symlog_scale(ax: Axes, *, x: bool = True, y: bool = True) -> None:
    """Put the requested axes on symlog scale for data spanning zero.

    Minor ticks are styled to be shorter than major ticks. If both x and y are
    scaled and matplotlib's automatic minor tick placement is inconsistent,
    both are disabled for uniformity.
    """
    axes_to_scale = []
    if x:
        ax.set_xscale("symlog")
        axes_to_scale.append(ax.xaxis)
    if y:
        ax.set_yscale("symlog")
        axes_to_scale.append(ax.yaxis)

    if not axes_to_scale:
        return

    _style_minor_ticks(ax)

    # If both axes are symlog-scaled, check for consistency
    if len(axes_to_scale) == 2:
        ax.figure.canvas.draw()

        x_minors = len(ax.xaxis.get_minor_locator()())
        y_minors = len(ax.yaxis.get_minor_locator()())

        if (x_minors > 0) != (y_minors > 0):
            ax.xaxis.set_minor_locator(plt.NullLocator())
            ax.yaxis.set_minor_locator(plt.NullLocator())


def apply_house_style() -> None:
    """Configure cnsplots with project house style: Cell palette, viridis density, Arial font."""
    # Remove Helvetica to prevent findfont warnings; Arial is installed
    cns.settings.font_sans_serif = ("Arial", "DejaVu Sans")
    cns.settings.panel_label_fontname = "Arial"

    # Fonttype 42 embeds TrueType, editable in Illustrator
    cns.settings.pdf_fonttype = 42

    cns.settings.panel_pad_left = PANEL_LABEL_PAD_LEFT_PX
    cns.settings.panel_pad_top = PANEL_LABEL_PAD_TOP_PX

    # Both are needed: setup_matplotlib sets the rcParams cycle for bare
    # matplotlib calls, palette_qual is what cns.* plot functions and
    # multipanel.panel() read when no explicit color_cycle is passed.
    cns.settings.palette_qual = HOUSE_PALETTE
    cns.settings.palette_seq = DENSITY_CMAP
    cns.setup_matplotlib(color_cycle=HOUSE_PALETTE)


def panel_labels(n: int) -> list[str]:
    """Return n unique panel labels: A-Z, then A1, A2, ... for n > 26."""
    return [
        chr(65 + i) if i < _ALPHABET_SIZE else f"A{i - _ALPHABET_SIZE + 1}"
        for i in range(n)
    ]


def figure_size_for_grid(
    n_rows: int, n_cols: int, shape: PanelShape = PanelShape.SQUARE
) -> tuple[int, int]:
    """Return the (width, height) page size in layout pixels that gives every cell ``shape``.

    The inverse of the usual arithmetic: the axes box is the fixed quantity and
    the page grows with the grid, so a 4-column figure holds panels the same
    size as a 1-column one instead of quartering them.

    An estimate only: the decoration constants are averages over real panels, and
    the actual reserve depends on how wide the tick labels render. ``fit_panels``
    corrects the remainder by measurement.
    """
    if n_rows < 1:
        raise ValueError(f"n_rows must be at least 1, got {n_rows}")
    if n_cols < 1:
        raise ValueError(f"n_cols must be at least 1, got {n_cols}")

    panel_width, panel_height = PANEL_SIZES_PX[shape]
    width = n_cols * (panel_width + PANEL_DECORATION_WIDTH_PX) + FIGURE_PAD_PX
    height = n_rows * (panel_height + PANEL_DECORATION_HEIGHT_PX) + FIGURE_PAD_PX

    if width > JOURNAL_WIDTH_PX:
        logger.warning(
            f"{n_cols} {shape} panels need {width} px, wider than the {JOURNAL_WIDTH_PX} px "
            f"journal page: the figure will be scaled down at typesetting. "
            f"Use fewer columns, or PanelShape.TALL for a narrower panel."
        )

    return width, height


def fit_panels(*, rect: tuple[float, float, float, float] | None = None) -> None:
    """Lay the current grid out and grow the page until each axes box matches its PanelShape.

    Replaces a bare ``tight_layout()`` at the end of a ``grid_axes`` figure. The
    estimate in ``figure_size_for_grid`` gets each panel within ~15 px, but the
    residual grows with the column count, because ``set_box_aspect`` shrinks a
    box whose cell is the wrong aspect and ``tight_layout`` does not give the
    space back. Text extents are only knowable after a draw, so the correction
    has to be measured rather than derived: each pass measures the shortfall,
    charges it across the grid, and re-lays out.

    ``rect`` is forwarded to ``tight_layout`` for callers reserving a strip
    (e.g. a figure-level legend).
    """
    fig = plt.gcf()
    spec: _GridSpecification | None = getattr(fig, _GRID_SPEC_ATTR, None)

    fig.tight_layout(rect=rect)
    if spec is None:
        return

    target_width, target_height = PANEL_SIZES_PX[spec.shape]
    points_per_px = 72 / fig.dpi

    for _ in range(_FIT_PASSES):
        fig.canvas.draw()
        visible = [ax for ax in spec.axes if ax.get_visible()]
        if not visible:
            return

        box = visible[0].get_window_extent()
        width_error = target_width - box.width * points_per_px
        height_error = target_height - box.height * points_per_px
        if abs(width_error) <= _FIT_TOLERANCE_PX and abs(height_error) <= _FIT_TOLERANCE_PX:
            return

        # One panel's shortfall is charged once per column (and per row), since
        # every cell in the grid is the same size and misses by the same amount.
        size = fig.get_size_inches()
        fig.set_size_inches(
            size[0] + spec.n_cols * width_error / 72,
            size[1] + spec.n_rows * height_error / 72,
        )
        fig.tight_layout(rect=rect)


def grid_axes(
    n_rows: int,
    n_cols: int,
    *,
    labels: Sequence[str] | None = None,
    shape: PanelShape = PanelShape.SQUARE,
    share_x: bool = False,
    share_y: bool = False,
) -> list[Axes]:
    """Create an n_rows x n_cols grid of aligned ``shape`` axes, row-major, and label each one.

    GridSpec allots equal cells regardless of how wide a panel's tick labels
    render, so column x and row y agree by construction -- which cns.multipanel
    cannot do, since it offsets every panel by its own measured decoration
    width. ``set_box_aspect`` then pins every axes box to the shape's aspect, so
    a grid whose row and column counts differ still holds equal panels.

    ``share_x``/``share_y`` put every panel on one common range and drop the
    interior tick labels, so a reader compares panels against one scale.

    Callers hide the cells they do not fill; ``tight_layout`` then reclaims the
    space. Run it after all panels are drawn, since it measures rendered text.
    """
    width_px, height_px = figure_size_for_grid(n_rows, n_cols, shape)
    panel_width, panel_height = PANEL_SIZES_PX[shape]

    cns.figure(width=width_px, height=height_px)
    fig = plt.gcf()
    grid = fig.subplots(n_rows, n_cols, squeeze=False, sharex=share_x, sharey=share_y)

    names = list(labels) if labels is not None else panel_labels(n_rows * n_cols)

    axes: list[Axes] = []
    for index in range(n_rows * n_cols):
        ax = grid[index // n_cols][index % n_cols]
        ax.set_box_aspect(panel_height / panel_width)
        if index < len(names):
            plt.sca(ax)
            cns.add_panel_label(names[index])
        axes.append(ax)

    setattr(
        fig,
        _GRID_SPEC_ATTR,
        _GridSpecification(n_rows=n_rows, n_cols=n_cols, shape=shape, axes=tuple(axes)),
    )

    return axes


def house_colors(indices: Sequence[int]) -> list[str]:
    """Return the house-palette hex colors at the given positions."""
    return cns.get_hexcolors_from_apalette(list(indices), HOUSE_PALETTE)


def observed_fitted_colors() -> tuple[str, str]:
    """Return the (observed, fitted) pair for measurement-versus-model figures."""
    observed, fitted = house_colors((_OBSERVED_INDEX, _FITTED_INDEX))
    return observed, fitted


def series_colors(n: int) -> list[str]:
    """Return n greyscale-separated house colors for overlaid series on one axes."""
    if n > len(_SERIES_INDICES):
        # Beyond the curated trio, luma separation is no longer guaranteed, so
        # fall back to palette order and let the caller carry the distinction on
        # a second channel (marker, dash) if print legibility matters.
        return house_colors(range(n))
    return house_colors(_SERIES_INDICES[:n])


def save_dual(stem: Path | str) -> None:
    """Save the current figure as a journal PDF plus a review PNG, given a stem with no extension."""
    # Append rather than replace suffixes: Path.with_suffix()/.stem would truncate
    # at the first dot, so a stem like "HD1328-4.YES0_corr" would silently collapse
    # to "HD1328-4" and two samples would overwrite each other's figures.
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)

    cns.savefig(stem.parent / f"{stem.name}.pdf")
    cns.savefig(stem.parent / f"{stem.name}.review.png")
