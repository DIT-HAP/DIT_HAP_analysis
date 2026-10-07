#!/usr/bin/env python3

"""
Per-View Fraction Distributions — Visualization
================================================

Three panels for one view, rendered from compute_view_fractions.py's two tables.
This is a pure renderer (ADR-0001): nothing is recomputed here, so the figure and
the tables can never disagree.

- A: per-term paralog fraction,
- B: groups per gene, with the view's cut marked — the quantile that defines moonlighting,
- C: per-term moonlighting fraction.

The two fraction panels pin x to [0, 1], a share's own range, so the same panel can
be read across the views; the breadth panel is in whole groups and takes the data's
own range, since a view with 1,400 terms reaches further than one with 170.

Input
-----
- --terms: view_fractions.tsv (paralog_fraction, moonlighting_fraction,
  moonlighting_cut_n_groups).
- --genes: view_gene_breadth.tsv (n_groups, is_moonlighting).
- --view: names the figure in the log line only; the file already lives in that
  view's folder.

Output
------
- --output: fraction_distributions.pdf (+ a .review.png sibling via save_dual).

Usage
-----
    python plot_fraction_distributions.py \\
        --terms results/3a_coherence/{dataset}/go_macrocomplex/view_fractions.tsv \\
        --genes results/3a_coherence/{dataset}/go_macrocomplex/view_gene_breadth.tsv \\
        --view go_macrocomplex \\
        --output results/3a_coherence/{dataset}/go_macrocomplex/fraction_distributions.pdf

Author:   Yusheng Yang (guidance) + Claude Sonnet 4.6 (implementation)
Date:     2026-10-05
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
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from coherence.fractions import MOONLIGHTING_QUANTILE  # noqa: E402
from figure_render.histogram import draw_histogram_panel  # noqa: E402
from figures import (  # noqa: E402
    FURNITURE_COLOR,
    PanelShape,
    apply_house_style,
    fit_panels,
    grid_axes,
    panel_labels,
    save_dual,
)
from io_table import read_file  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Bins per panel. The fraction panels are shares on a pinned 0-1 range, so their
# bin width is fixed across views (0.05); the breadth panel takes 30 equal bins
# over the view's observed range instead — one bin per integer would be a comb at
# the 100 px panel width, and the tail runs to ~120 groups.
_FRACTION_BINS = 20
_BREADTH_BINS = 30

_CUT_LINEWIDTH = 1.0


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class DistributionPlotConfig:
    """The view's two computed tables plus the figure path."""
    terms: Path
    genes: Path
    view: str
    output: Path

    def validate(self) -> None:
        """Raise ValueError if an input is missing, then create the output dir."""
        for path in [self.terms, self.genes]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
def plot_distributions(terms: pd.DataFrame, genes: pd.DataFrame) -> None:
    """Draw the three panels of one view."""
    apply_house_style()
    axes = grid_axes(1, 3, labels=panel_labels(3), shape=PanelShape.SQUARE)

    draw_histogram_panel(
        axes[0], terms["paralog_fraction"], bins=_FRACTION_BINS,
        xlabel="Paralog fraction", ylabel="# terms", title="Paralog buffering",
    )
    # The cut the moonlighting column is taken with. Drawn as furniture, not as a
    # fitted quantity: it is the quantile of the panel's own data.
    cut = float(terms["moonlighting_cut_n_groups"].iloc[0])
    draw_histogram_panel(
        axes[1], genes["n_groups"], bins=_BREADTH_BINS,
        xlabel="Groups per gene", ylabel="# genes", title="Gene breadth",
    )
    axes[1].axvline(cut, color=FURNITURE_COLOR, linestyle="--", linewidth=_CUT_LINEWIDTH,
                    label=f"P{MOONLIGHTING_QUANTILE * 100:g} = {cut:g}")
    axes[1].legend(loc="upper right")

    draw_histogram_panel(
        axes[2], terms["moonlighting_fraction"], bins=_FRACTION_BINS,
        xlabel="Moonlighting fraction", ylabel="# terms", title="Moonlighting",
    )
    # A share's own range, pinned so the same panel can be read across views.
    for ax in (axes[0], axes[2]):
        ax.set_xlim(0.0, 1.0)

    fit_panels()


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: DistributionPlotConfig) -> None:
    """Read the computed tables, draw the three distributions, save the figure."""
    config.validate()
    terms = read_file(config.terms)
    genes = read_file(config.genes)
    if terms.empty or genes.empty:
        logger.warning(f"[{config.view}] empty fractions table; no figure written")
        return

    plot_distributions(terms, genes)
    save_dual(config.output.with_suffix(""))
    logger.success(
        f"[{config.view}] {len(terms):,} terms / {len(genes):,} genes; "
        f"cut = {float(terms['moonlighting_cut_n_groups'].iloc[0]):g} groups/gene; wrote {config.output}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Plot one view's paralog / moonlighting / breadth distributions")
    parser.add_argument("--terms", type=Path, required=True, help="view_fractions.tsv")
    parser.add_argument("--genes", type=Path, required=True, help="view_gene_breadth.tsv")
    parser.add_argument("--view", required=True, help="The view's name, for the log line")
    parser.add_argument("--output", type=Path, required=True, help="Output PDF (a .review.png sibling is written too)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render the figure, report the path."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = DistributionPlotConfig(
            terms=args.terms, genes=args.genes, view=args.view, output=args.output,
        )
        run(config)
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
