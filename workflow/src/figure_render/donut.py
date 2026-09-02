"""Single-axes donut primitive built on ``cns.donutplot``.

``cns.donutplot`` hard-codes the column name as both the center annotation and
the legend title, which is meaningless for part/whole donuts driven by two raw
counts. This module owns the post-patching: center text replaced by caller
supplied text and the legend title dropped.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-28
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import cnsplots as cns
import matplotlib.axes
import matplotlib.legend
import matplotlib.text
import pandas as pd

# =============================================================================
# CONSTANTS
# =============================================================================
STATUS_COLUMN = "status"

# Gap between the ring's lowest point and the legend's top, in layout pixels.
# cns.donutplot anchors a bottom legend at (0.5, -0.05) in axes fraction, but
# the ring only spans |y| <= 1.0 inside ylim +/-1.25, so ~10% of the axes height
# is already empty below the ring before that offset applies — ~20 layout px on
# a 110 px panel, which reads as a legend detached from its own donut.
LEGEND_RING_GAP_PX = 7

# Axes fraction of the ring's lowest point: radius 1.0 within ylim +/-1.25.
_RING_BOTTOM_FRACTION = 0.1

# The anchor places the legend's padded box, whose top sits this far above the
# rendered handle/label extent. Measured at 3.5 layout px for the house style
# (borderpad 0.4 at 7 pt); subtracted so LEGEND_RING_GAP_PX means the visible gap.
_LEGEND_BORDER_PAD_PX = 3.5


# =============================================================================
# CORE LOGIC
# =============================================================================
def draw_donut_panel(
    ax: matplotlib.axes.Axes,
    part: int,
    whole: int,
    *,
    part_label: str,
    whole_label: str,
    center_text: str,
) -> matplotlib.text.Text:
    """Draw one part/whole donut on ``ax`` with ``center_text`` in the hole.

    Returns the center annotation so callers can style or assert on it. The
    ring starts at 12 o'clock and runs clockwise with ``part_label``'s wedge
    first, legend beneath. Counts of zero on both sides are handled by the
    caller — this primitive always draws a wedge ring, however degenerate.
    """
    counts = [part, whole]
    status_df = pd.DataFrame(
        {STATUS_COLUMN: [part_label] * part + [whole_label] * whole}
    )
    cns.donutplot(
        data=status_df,
        x=STATUS_COLUMN,
        legend="bottom",
        order=[part_label, whole_label],
        ax=ax,
    )

    # cns.donutplot puts the column name ("status") at the center and as the
    # legend title; replace the former with the caller's text and drop the
    # latter, which would otherwise repeat the meaningless column name. The
    # panel-label Text may precede the center annotation in ax.texts, so pick
    # the annotation by its original column-name content instead of position.
    center = next(text for text in ax.texts if text.get_text() == STATUS_COLUMN)
    center.set_text(center_text)
    center.set_fontsize(cns.settings.legend_fontsize)

    legend = ax.get_legend()
    legend.set_title(None)
    _tighten_bottom_legend(ax, legend)
    return center


def _tighten_bottom_legend(
    ax: matplotlib.axes.Axes, legend: matplotlib.legend.Legend
) -> None:
    """Re-anchor a bottom legend LEGEND_RING_GAP_PX below the ring, not below the axes box."""
    axes_height_px = ax.get_position().height * ax.figure.get_size_inches()[1] * 72
    if axes_height_px <= 0:
        return

    offset_px = LEGEND_RING_GAP_PX - _LEGEND_BORDER_PAD_PX
    legend.set_bbox_to_anchor(
        (0.5, _RING_BOTTOM_FRACTION - offset_px / axes_height_px),
        transform=ax.transAxes,
    )
