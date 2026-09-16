"""
Spike-In Panel Renderer
=======================

Panel (d) of the PCR QC figure: spike-in dilution linearity. Thin wrapper over
``cns.regplot`` — the scatter, the linear fit, the Pearson r/P annotation and
the legend are all cnsplots' own; this module only excludes the zero-dilution
reference and fixes the axis labels.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-02
Version:  3.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import matplotlib.pyplot as plt
import pandas as pd
import cnsplots as cns
from loguru import logger

# The zero-dilution reference sample: no spike-in reads to be linear against.
REFERENCE_SAMPLE = "Spikein0"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def render_spikein_panel(
    ax: plt.Axes,
    spikein: pd.DataFrame,
    *,
    marker_size: float = 20,
    add_equation: bool = True,
    hue_order: list[str] | None = None,
) -> None:
    """Render the spike-in dilution linearity panel via ``cns.regplot``.

    Draws log2(relative dilution ratio) against log2(relative read ratio), one
    colour per spike-in insertion site, with a single linear fit and Pearson
    correlation over all points.

    Parameters
    ----------
    ax : plt.Axes
        Target axes, already styled by ``figures.apply_house_style``.
    spikein : pd.DataFrame
        Spike-in stats with columns Sample, Name, Relative_Dilution_Ratio,
        Relative_Read_Ratio.
    marker_size : float, default 20
        Scatter marker size in points², forwarded to ``cns.regplot``'s ``s``.
    add_equation : bool, default True
        Annotate the fitted equation and R² in the bottom-right.
    hue_order : list of str, optional
        Explicit ``Name`` order, so the site→colour mapping does not depend on
        pandas' group ordering. Defaults to the encountered order.

    Notes
    -----
    ``Spikein0`` is dropped: it is the zero-dilution reference, so its reads say
    nothing about linearity, and its floored read count would anchor the fit at
    -inf. Matches the source notebook.
    """
    spikein_filtered = spikein.query("Sample != @REFERENCE_SAMPLE").copy()
    if spikein_filtered.empty:
        raise ValueError(f"No spike-in rows left after dropping {REFERENCE_SAMPLE}")

    if hue_order is None:
        hue_order = sorted(spikein_filtered["Name"].unique())

    # color=<column> (not hue=): one legend entry per insertion site, but a
    # single overall fit — hue= would fit and annotate five separate lines.
    cns.regplot(
        data=spikein_filtered,
        x="Relative_Dilution_Ratio",
        y="Relative_Read_Ratio",
        color="Name",
        hue_order=hue_order,
        s=marker_size,
        add_equation=add_equation,
        ax=ax,
    )

    ax.set_xlabel(r"$\log_2$(relative dilution ratio)")
    ax.set_ylabel(r"$\log_2$(relative read ratio)")

    # cnsplots moves the legend to the right margin once add_equation is on —
    # the only free space, since the series fill the diagonal. Re-run it with an
    # empty title so the colour column's name does not leak in as "Name".
    if ax.get_legend() is not None:
        cns.take_legend_out(title="", ax=ax)

    logger.debug(
        f"Rendered spike-in panel: {len(spikein_filtered)} points, "
        f"{len(hue_order)} insertion sites"
    )
