"""Colour and ordering vocabulary for the coherence stage's cross-source figures.

Two plotting scripts now draw more than one source at once — the coherence
comparison figure (all sources + the de-duplicated representative set) and the
pooled attribution figure's label-frequency bars. This module exists so they
cannot drift apart on what colour `go_bp` is.

It deliberately does NOT live in `sources.py`: that module is imported by the
compute stage, and a colour helper there would drag matplotlib into a
compute-only environment. The palette lookup is likewise deferred to call time,
because `house_colors` reads whatever palette `apply_house_style()` installed.
"""
# =============================================================================
# IMPORTS
# =============================================================================
from __future__ import annotations

from collections.abc import Sequence


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Draw/legend order. The five real sources follow config.coherence.sources, so the
# figures read in the same order as the pipeline's own registry; `dedup` is the
# de-duplicated representative subset, which only the comparison figure draws.
SOURCE_ORDER = [
    "go_macrocomplex", "go_cc", "go_bp", "kegg_brite", "kegg_pathway", "dedup",
]

# Cell palette positions, hand-picked rather than range(6): adjacent palette slots
# 0/1 (red, teal) differ by 0.003 in relative luminance and would print as a single
# tone. The first three match plot_redundancy_network.py's assignment, so a source
# keeps the same colour between the network page and these figures.
_SOURCE_PALETTE_INDICES = {
    "go_bp": 2, "go_cc": 1, "go_macrocomplex": 0,
    "kegg_brite": 3, "kegg_pathway": 6, "dedup": 4,
}


# =============================================================================
# CORE LOGIC
# =============================================================================
def source_colors(sources: Sequence[str] | None = None) -> dict[str, str]:
    """Return {source: house-palette hex} for `sources` (default: SOURCE_ORDER)."""
    # A source the map does not know takes the palette slot at its own position;
    # config, not this file, decides which sources exist, so an unmapped one should
    # still draw (in a distinct colour) rather than raise.
    from figures import house_colors  # noqa: PLC0415 - see the module docstring

    labels = list(sources) if sources is not None else list(SOURCE_ORDER)
    indices = [
        _SOURCE_PALETTE_INDICES.get(label, position)
        for position, label in enumerate(labels)
    ]
    return dict(zip(labels, house_colors(indices)))
