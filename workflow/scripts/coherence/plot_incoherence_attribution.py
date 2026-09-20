#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coherence Incoherence Attribution — Visualization
=================================================

Renders the figure for compute_incoherence_attribution.py: a grid of DR-DL
scatter panels, one per top-N incoherent group, with members coloured by GMM
component, plus a closing panel counting the attribution labels.

This is the plotting companion, split out per ADR-0001. It re-fits nothing: the
GMM components and the per-member coordinates it draws were computed once and
persisted to incoherence_split_points.parquet, so the figure is a pure renderer
and the components it shows are the same ones the attribution table was labelled
from.

Input
-----
- incoherence_attribution.tsv: per-group attribution rows. Only is_incoherent
  rows are drawn, and the label counts come from all of them.
- incoherence_split_points.parquet: one row per member of every incoherent group
  (group_id, Systematic ID, norm_DR, norm_DL, component).

Output
------
- incoherence_attribution.pdf (+ a .review.png sibling via save_dual)

Usage
-----
    python plot_incoherence_attribution.py \\
        --table results/3a_coherence/{dataset}/go_macrocomplex/incoherence_attribution.tsv \\
        --points results/3a_coherence/{dataset}/go_macrocomplex/incoherence_split_points.parquet \\
        --top-n-plot 16 \\
        --output results/3a_coherence/{dataset}/go_macrocomplex/incoherence_attribution.pdf

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-23
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import math
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import cnsplots as cns
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from figures import (  # noqa: E402
    FURNITURE_COLOR,
    PanelShape,
    apply_house_style,
    fit_panels,
    grid_axes,
    house_colors,
    panel_labels,
    save_dual,
)
from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Fixed category order, used as hue_order on every panel: without it seaborn
# assigns colours in order of first appearance, so "core" would take a different
# colour in a panel that happens to lack "minor".
_COMPONENT_ORDER = ["core", "minor", "single"]

# Scatter panels per row. 4 keeps a 16-group figure to 5 rows; the resulting page
# is wider than the one-page cap (grid_axes warns), which is accepted for what is
# a working diagnostic figure rather than a typeset one.
_MAX_COLUMNS = 4

# Group names are shortened to this many characters for a panel title: a 100 px
# SQUARE panel fits about 25 at the house title size, minus the ellipsis.
_TITLE_NAME_WIDTH = 22

_REQUIRED_TABLE_COLUMNS = ["group_id", "group_name", "median_pairwise_distance_z",
                           "is_incoherent", "attribution_label"]
_REQUIRED_POINT_COLUMNS = ["group_id", "norm_DR", "norm_DL", "component"]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Inputs, outputs, and parameters for the attribution figure."""
    table: Path
    points: Path
    output: Path
    top_n_plot: int = 16

    def validate(self) -> None:
        """Raise ValueError if inputs are missing or params invalid, then make output dirs."""
        for path in [self.table, self.points]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.top_n_plot < 1:
            raise ValueError(f"top_n_plot must be >= 1: {self.top_n_plot}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
def component_colors() -> dict[str, str]:
    """Marker colour per GMM component, resolved from the house palette at call time."""
    # Palette positions 1 and 2 (teal, amber) are the pair with real greyscale
    # separation — relative luminances 0.18 vs 0.42. The adjacent positions 0 and 1
    # (red, teal) sit at 0.178 and 0.175, i.e. they would print as the same tone, and
    # core-vs-minor is the one distinction these panels exist to show.
    core, minor = house_colors((1, 2))
    return {"core": core, "minor": minor, "single": FURNITURE_COLOR}


def _panel_title(row: pd.Series) -> str:
    """Panel title: group name, its z-score, and the GMM silhouette when there was one."""
    # Three lines, each short enough for a SQUARE panel. A 100 px panel holds about
    # 25 characters at the 8 pt house title size, so a two-line
    # "name / z=3.33, sil=0.76 [conditional_module]" overflows into the neighbouring
    # panel — measured on the real go_macrocomplex figure, where line 2 is 37
    # characters. The name is shortened on word boundaries.
    name = textwrap.shorten(str(row["group_name"]), width=_TITLE_NAME_WIDTH, placeholder="…")
    silhouette = row["gmm_silhouette"] if "gmm_silhouette" in row else float("nan")
    stats = f"z={row['median_pairwise_distance_z']:.2f}"
    if pd.notna(silhouette):
        stats += f", sil={silhouette:.2f}"
    return f"{name}\n{stats}\n[{row['attribution_label']}]"


def plot_attribution(table: pd.DataFrame, points: pd.DataFrame, top_n: int) -> None:
    """Draw the attribution figure onto a fresh house-styled grid."""
    # The leading square panel and the square `.review.png` sibling come from
    # `save_dual`, called by the caller. No incoherent groups renders a single
    # placeholder panel rather than an empty figure.
    apply_house_style()

    incoherent = table[table["is_incoherent"]].head(top_n)
    n_panels = len(incoherent) + 1  # + the label-frequency panel
    n_cols = min(_MAX_COLUMNS, n_panels)
    n_rows = math.ceil(n_panels / n_cols)
    axes = grid_axes(n_rows, n_cols, labels=panel_labels(n_panels), shape=PanelShape.SQUARE)
    for ax in axes[n_panels:]:
        # grid_axes fills every cell; fit_panels measures only visible axes, so
        # the unfilled ones must be hidden (not just deleted) to be discounted.
        ax.set_visible(False)

    if incoherent.empty:
        axes[0].text(
            0.5, 0.5, "No incoherent groups (z > threshold)",
            ha="center", va="center", transform=axes[0].transAxes,
        )
        axes[0].set_axis_off()
        fit_panels()
        return

    colors = component_colors()
    palette = [colors[component] for component in _COMPONENT_ORDER]

    for index, (_, row) in enumerate(incoherent.iterrows()):
        group_points = points[points["group_id"] == row["group_id"]]
        ax = axes[index]
        if group_points.empty:
            ax.text(0.5, 0.5, "No fitted members", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
        else:
            cns.scatterplot(
                group_points, "norm_DR", "norm_DL",
                hue="component", hue_order=_COMPONENT_ORDER,
                palette=palette, legend=(index == 0), ax=ax,
            )
            ax.set_xlabel("norm DR")
            ax.set_ylabel("norm DL/10")
        ax.set_title(_panel_title(row))

    # Label-frequency panel: a plain count bar. cns.barplot aggregates a mean per
    # category with optional significance testing, which is not this.
    counts = table[table["is_incoherent"]]["attribution_label"].value_counts().sort_values()
    ax_freq = axes[len(incoherent)]
    ax_freq.barh(counts.index, counts.values, color=house_colors((3,))[0])
    ax_freq.set_xlabel("Number of incoherent groups")
    ax_freq.set_title("Label frequency")

    fit_panels()


@logger.catch(reraise=True)
def run(config: PlotConfig) -> None:
    """Load -> draw -> save the attribution figure."""
    config.validate()
    table = pd.read_csv(config.table, sep="\t")
    missing = [col for col in _REQUIRED_TABLE_COLUMNS if col not in table.columns]
    if missing:
        raise ValueError(f"attribution table missing required column(s) {missing} (have: {list(table.columns)})")
    points = read_parquet(config.points)
    if not points.empty:
        missing = [col for col in _REQUIRED_POINT_COLUMNS if col not in points.columns]
        if missing:
            raise ValueError(f"split points missing required column(s) {missing} (have: {list(points.columns)})")

    n_incoherent = int(table["is_incoherent"].sum()) if not table.empty else 0
    logger.info(
        f"{len(table):,} groups, {n_incoherent:,} incoherent, {len(points):,} split points; "
        f"drawing top {min(config.top_n_plot, n_incoherent)}"
    )

    plot_attribution(table, points, config.top_n_plot)
    save_dual(config.output.with_suffix(""))
    logger.success(f"Wrote {config.output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Plot the incoherence attribution figure")
    parser.add_argument("--table", type=Path, required=True, help="incoherence_attribution.tsv")
    parser.add_argument("--points", type=Path, required=True, help="incoherence_split_points.parquet")
    parser.add_argument("--top-n-plot", type=int, default=16, help="How many top-incoherent groups to scatter")
    parser.add_argument("--output", type=Path, required=True, help="Output attribution figure PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run plotting, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotConfig(
            table=args.table,
            points=args.points,
            output=args.output,
            top_n_plot=args.top_n_plot,
        )
        run(config)
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
