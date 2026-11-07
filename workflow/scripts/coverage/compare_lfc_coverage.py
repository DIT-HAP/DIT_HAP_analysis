#!/usr/bin/env python3

"""
Compare LFC Coverage
====================

Cross-dataset half of the LFC-coverage comparison: stack the per-dataset
lfc_coverage_stats.tsv tables in config order, write the combined TSV, then
draw the grouped bar charts from that same table (so a figure cannot disagree
with the table it ships beside).

Three figures share the data:

- ``lfc_coverage_by_fypo`` — vertical bars, FYPO class on x, all datasets.
- ``lfc_coverage_by_fypo_horizontal`` — the same chart with x and y swapped
  (FYPO class down the y axis, coverage across x).
- ``lfc_coverage_by_fypo_two_libraries`` — the swapped layout restricted to
  LD_haploid and HD_DIT_HAP, displayed as "Traditional Tn-Seq" and "DIT_HAP".

Every chart labels each bar with its absolute LFC-covered gene count. The
three-dataset charts star the `total` group's LD_DIT_HAP bar and footnote the
exact counts, because LD_DIT_HAP covers fewer genes overall than LD_haploid —
the essential-coverage advantage does not extend to the total. No confidence
interval is drawn: the LFC table is already the output of the hard filter +
biological-replicate merge, so a binomial interval over genes would not describe
replicate uncertainty.

Input
-----
- results/2a_coverage/{dataset}/lfc_coverage_stats.tsv (one per dataset)

Output
------
- results/2a_coverage/cross_dataset/lfc_coverage_comparison.tsv
- results/2a_coverage/cross_dataset/lfc_coverage_by_fypo.{pdf,review.png}
- results/2a_coverage/cross_dataset/lfc_coverage_by_fypo_horizontal.{pdf,review.png}
- results/2a_coverage/cross_dataset/lfc_coverage_by_fypo_two_libraries.{pdf,review.png}

Usage
-----
    python workflow/scripts/coverage/compare_lfc_coverage.py \\
        --stats results/2a_coverage/LD_haploid/lfc_coverage_stats.tsv \\
                results/2a_coverage/HD_DIT_HAP/lfc_coverage_stats.tsv \\
                results/2a_coverage/LD_DIT_HAP/lfc_coverage_stats.tsv \\
        --output-tsv results/2a_coverage/cross_dataset/lfc_coverage_comparison.tsv \\
        --output-figure results/2a_coverage/cross_dataset/lfc_coverage_by_fypo.pdf \\
        --output-figure-horizontal \\
        results/2a_coverage/cross_dataset/lfc_coverage_by_fypo_horizontal.pdf \\
        --output-figure-two-libraries \\
        results/2a_coverage/cross_dataset/lfc_coverage_by_fypo_two_libraries.pdf

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-07
Version:  1.1.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
import sys

# 2. Third-party Imports
import cnsplots as cns
from loguru import logger
from matplotlib.axes import Axes
import numpy as np
import pandas as pd

# 3. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from coverage.core import LFC_TOTAL_LABEL, combine_lfc_coverage_stats
from figures import (
    FIT_ANNOTATION_FONT_PT,
    FURNITURE_COLOR,
    apply_house_style,
    house_colors,
    save_dual,
    series_colors,
)
from logging_setup import setup_logger

# =============================================================================
# CONSTANTS
# =============================================================================
REQUIRED_STATS_COLUMNS = (
    "dataset",
    "fypo_class",
    "n_total",
    "n_dr",
    "n_lfc",
    "pct_dr",
    "pct_lfc",
)

# Display labels for the axis; the stats table keeps the raw class names.
CLASS_LABELS = {
    "inviable": "Inviable",
    "viable": "Viable",
    "depends_on_conditions": "Depends on\nconditions",
    "unknown": "Unknown",
    LFC_TOTAL_LABEL: "Total",
}

PERCENT_PER_UNIT = 100.0
BAR_GROUP_WIDTH = 0.8
# Headroom above 100% for the count labels; the plotted data is still the
# untruncated 0-100% range, and the ticks stop at 100.
AXIS_HEADROOM_PCT = 108
X_AXIS_LABEL = "FYPO viability class"
Y_AXIS_LABEL = "LFC-covered genes (%)"
TITLE = "Gene coverage by FYPO viability class"

# The total row's caveat: LD_DIT_HAP's all-gene coverage sits below LD_haploid's,
# so the aggregate bar is starred and the exact counts are footnoted. Both only
# appear in the charts that actually plot those two datasets.
STAR_DATASET = "LD_DIT_HAP"
REFERENCE_DATASET = "LD_haploid"
STAR_CLASS = LFC_TOTAL_LABEL

# The two-library chart's display names. The stats tables keep the registry
# names; only the figure legend is renamed.
TWO_LIBRARY_DATASETS = ("LD_haploid", "HD_DIT_HAP")
TWO_LIBRARY_LABELS = {"LD_haploid": "Traditional Tn-Seq", "HD_DIT_HAP": "DIT_HAP"}

# The essential class is the figure's message; mark it with the palette's first
# (accent) colour as a background band, without changing the dataset colour code.
FOCUS_CLASS = "inviable"
FOCUS_BAND_ALPHA = 0.08


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class FigureLayout:
    """Panel/canvas geometry plus the axes-fraction legend and footnote positions."""

    panel_width_pt: int
    panel_height_pt: int
    canvas_width_pt: int
    margin_bottom_pt: int
    legend_y: float
    footnote_y: float


# Fifteen bars plus a 4-digit count label per bar need more width than the
# standard 150 pt WIDE panel; 230 pt keeps the three bars of a group legible.
# The horizontal layouts are taller (five class groups down the y axis) and need
# a wider canvas for the y tick labels ("Traditional Tn-Seq" is the longest).
VERTICAL_LAYOUT = FigureLayout(
    panel_width_pt=230,
    panel_height_pt=110,
    canvas_width_pt=310,
    margin_bottom_pt=70,
    legend_y=-0.40,
    footnote_y=-0.56,
)
HORIZONTAL_LAYOUT = FigureLayout(
    panel_width_pt=230,
    panel_height_pt=145,
    canvas_width_pt=380,
    margin_bottom_pt=55,
    legend_y=-0.20,
    footnote_y=-0.32,
)


@dataclass(kw_only=True, slots=True, frozen=True)
class CompareLfcCoverageConfig:
    """Per-dataset stats tables, the combined TSV, and the three figure paths."""

    stats: list[Path]
    output_tsv: Path
    output_figure: Path
    output_figure_horizontal: Path | None = None
    output_figure_two_libraries: Path | None = None

    def validate(self) -> None:
        """Raise ValueError if an input is missing, then ensure output dirs exist."""
        for path in self.stats:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for output in [
            self.output_tsv,
            self.output_figure,
            self.output_figure_horizontal,
            self.output_figure_two_libraries,
        ]:
            if output is not None:
                output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
def read_stats_tables(paths: list[Path]) -> pd.DataFrame:
    """Read and stack the per-dataset stats tables in command-line order."""
    tables = [pd.read_csv(path, sep="\t") for path in paths]
    for path, table in zip(paths, tables, strict=True):
        missing_columns = set(REQUIRED_STATS_COLUMNS) - set(table.columns)
        if missing_columns:
            raise ValueError(f"{path} is missing required columns: {sorted(missing_columns)}")
    return combine_lfc_coverage_stats(tables)


def _pivot(
    comparison: pd.DataFrame, column: str, classes: list[str], datasets: Sequence[str]
) -> pd.DataFrame:
    """Pivot one stats column to a FYPO-class x dataset frame in display order."""
    # aggfunc="first": (fypo_class, dataset) is unique by contract, so this is a
    # reshape rather than an aggregation, and it keeps n_lfc an integer column.
    return comparison.pivot_table(
        index="fypo_class", columns="dataset", values=column, aggfunc="first"
    ).loc[classes, list(datasets)]


def _total_lfc_count(comparison: pd.DataFrame, dataset: str) -> int:
    """LFC-covered gene count for one dataset's total row."""
    row = comparison[
        (comparison["dataset"] == dataset) & (comparison["fypo_class"] == LFC_TOTAL_LABEL)
    ]
    return int(row["n_lfc"].iloc[0])


def _bar_geometry(classes: list[str], datasets: Sequence[str]) -> tuple[np.ndarray, float]:
    """Group centres and per-series bar width shared by both orientations."""
    positions = np.arange(len(classes))
    bar_width = BAR_GROUP_WIDTH / len(datasets)
    group_offsets = positions - BAR_GROUP_WIDTH / 2 + bar_width / 2
    return group_offsets, bar_width


def _count_label(dataset: str, fypo_class: str, count: float) -> str:
    """Bar label: the absolute LFC-covered count, starred on the caveat bar."""
    suffix = "*" if (dataset == STAR_DATASET and fypo_class == STAR_CLASS) else ""
    return f"{int(count)}{suffix}"


def _draw_vertical_bars(
    ax: Axes,
    pct: pd.DataFrame,
    counts: pd.DataFrame,
    datasets: Sequence[str],
    classes: list[str],
    colors: list[str],
    labels: Mapping[str, str],
) -> None:
    """One vertical bar + rotated count label per dataset within every FYPO group."""
    group_offsets, bar_width = _bar_geometry(classes, datasets)

    for series_index, dataset in enumerate(datasets):
        offsets = group_offsets + series_index * bar_width
        values = pct[dataset].to_numpy()
        ax.bar(
            offsets,
            values,
            width=bar_width,
            color=colors[series_index],
            linewidth=0,
            label=labels.get(dataset, dataset),
        )
        for offset, value, count, fypo_class in zip(
            offsets, values, counts[dataset].to_numpy(), classes, strict=True
        ):
            ax.text(
                offset,
                value + 1.0,
                _count_label(dataset, fypo_class, count),
                rotation=90,
                ha="center",
                va="bottom",
                fontsize=FIT_ANNOTATION_FONT_PT,
                color=FURNITURE_COLOR,
            )


def _draw_horizontal_bars(
    ax: Axes,
    pct: pd.DataFrame,
    counts: pd.DataFrame,
    datasets: Sequence[str],
    classes: list[str],
    colors: list[str],
    labels: Mapping[str, str],
) -> None:
    """One horizontal bar + count label per dataset within every FYPO group."""
    group_offsets, bar_width = _bar_geometry(classes, datasets)

    for series_index, dataset in enumerate(datasets):
        offsets = group_offsets + series_index * bar_width
        values = pct[dataset].to_numpy()
        ax.barh(
            offsets,
            values,
            height=bar_width,
            color=colors[series_index],
            linewidth=0,
            label=labels.get(dataset, dataset),
        )
        for offset, value, count, fypo_class in zip(
            offsets, values, counts[dataset].to_numpy(), classes, strict=True
        ):
            ax.text(
                value + 1.0,
                offset,
                _count_label(dataset, fypo_class, count),
                rotation=0,
                ha="left",
                va="center",
                fontsize=FIT_ANNOTATION_FONT_PT,
                color=FURNITURE_COLOR,
            )


def _focus_class(ax: Axes, classes: list[str], accent: str, *, horizontal: bool) -> None:
    """Shade the essential class and colour its tick label with the accent."""
    focus_index = classes.index(FOCUS_CLASS)
    span = (focus_index - 0.5, focus_index + 0.5)
    if horizontal:
        ax.axhspan(*span, color=accent, alpha=FOCUS_BAND_ALPHA, linewidth=0, zorder=0)
        tick = ax.get_yticklabels()[focus_index]
    else:
        ax.axvspan(*span, color=accent, alpha=FOCUS_BAND_ALPHA, linewidth=0, zorder=0)
        tick = ax.get_xticklabels()[focus_index]
    tick.set_color(accent)
    tick.set_fontweight("bold")


def render_comparison_figure(
    comparison: pd.DataFrame,
    output_stem: Path,
    *,
    layout: FigureLayout,
    horizontal: bool = False,
    datasets: Sequence[str] | None = None,
    dataset_labels: Mapping[str, str] | None = None,
) -> None:
    """Draw one grouped LFC-coverage figure in the requested orientation."""
    apply_house_style()

    available = list(dict.fromkeys(comparison["dataset"]))
    selected = list(datasets) if datasets is not None else available
    classes = list(dict.fromkeys(comparison["fypo_class"]))
    labels = dataset_labels or {}
    colors = series_colors(len(selected))
    accent = house_colors([0])[0]

    pct = _pivot(comparison, "pct_lfc", classes, selected) * PERCENT_PER_UNIT
    counts = _pivot(comparison, "n_lfc", classes, selected)

    multipanel = cns.multipanel(max_width=layout.canvas_width_pt, unit="pt")
    ax = multipanel.panel(
        label="",
        width=layout.panel_width_pt,
        height=layout.panel_height_pt,
        margin_bottom=layout.margin_bottom_pt,
    )

    if horizontal:
        _draw_horizontal_bars(ax, pct, counts, selected, classes, colors, labels)
        ax.set_yticks(np.arange(len(classes)), [CLASS_LABELS[c] for c in classes])
        ax.invert_yaxis()
        # After set_yticks: it recreates the tick Text objects, so styling them
        # before would be discarded. The band's zorder keeps it behind the bars.
        _focus_class(ax, classes, accent, horizontal=True)
        ax.set_xlabel(Y_AXIS_LABEL)
        ax.set_ylabel(X_AXIS_LABEL)
        ax.set_xlim(0, AXIS_HEADROOM_PCT)
        ax.set_xticks([0, 20, 40, 60, 80, 100])
    else:
        _draw_vertical_bars(ax, pct, counts, selected, classes, colors, labels)
        ax.set_xticks(np.arange(len(classes)), [CLASS_LABELS[c] for c in classes])
        _focus_class(ax, classes, accent, horizontal=False)
        ax.set_xlabel(X_AXIS_LABEL)
        ax.set_ylabel(Y_AXIS_LABEL)
        ax.set_ylim(0, AXIS_HEADROOM_PCT)
        ax.set_yticks([0, 20, 40, 60, 80, 100])

    ax.set_title(TITLE)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, layout.legend_y),
        ncol=len(selected),
        frameon=False,
    )

    if STAR_DATASET in selected and REFERENCE_DATASET in selected:
        footnote = (
            f"{STAR_DATASET} {_total_lfc_count(comparison, STAR_DATASET):,} < "
            f"{REFERENCE_DATASET} {_total_lfc_count(comparison, REFERENCE_DATASET):,}"
        )
        ax.text(
            0.5,
            layout.footnote_y,
            footnote,
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=FIT_ANNOTATION_FONT_PT,
            color=FURNITURE_COLOR,
        )

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)


def run(config: CompareLfcCoverageConfig) -> None:
    """Stack the stats tables, write the comparison TSV, then draw the figures."""
    config.validate()

    comparison = read_stats_tables(config.stats)
    comparison.to_csv(config.output_tsv, sep="\t", index=False)
    logger.info(
        f"Wrote {len(comparison)} rows across "
        f"{comparison['dataset'].nunique()} datasets to {config.output_tsv}"
    )

    render_comparison_figure(
        comparison,
        config.output_figure.with_suffix(""),
        layout=VERTICAL_LAYOUT,
    )
    if config.output_figure_horizontal is not None:
        render_comparison_figure(
            comparison,
            config.output_figure_horizontal.with_suffix(""),
            layout=HORIZONTAL_LAYOUT,
            horizontal=True,
        )
    if config.output_figure_two_libraries is not None:
        render_comparison_figure(
            comparison,
            config.output_figure_two_libraries.with_suffix(""),
            layout=HORIZONTAL_LAYOUT,
            horizontal=True,
            datasets=TWO_LIBRARY_DATASETS,
            dataset_labels=TWO_LIBRARY_LABELS,
        )

    logger.success("LFC coverage figures complete")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Compare LFC coverage across datasets")
    parser.add_argument(
        "--stats",
        type=Path,
        nargs="+",
        required=True,
        help="Per-dataset lfc_coverage_stats.tsv files, in display order",
    )
    parser.add_argument(
        "--output-tsv", type=Path, required=True, help="Output combined comparison TSV"
    )
    parser.add_argument("--output-figure", type=Path, required=True, help="Output figure PDF path")
    parser.add_argument(
        "--output-figure-horizontal",
        type=Path,
        help="Optional horizontal (swapped-axis) figure PDF path",
    )
    parser.add_argument(
        "--output-figure-two-libraries",
        type=Path,
        help="Optional two-library figure PDF path (Traditional Tn-Seq vs DIT_HAP)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging"
    )
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, write the table, render the figures."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = CompareLfcCoverageConfig(
            stats=args.stats,
            output_tsv=args.output_tsv,
            output_figure=args.output_figure,
            output_figure_horizontal=args.output_figure_horizontal,
            output_figure_two_libraries=args.output_figure_two_libraries,
        )
        run(config)
    except Exception as e:
        logger.exception(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
