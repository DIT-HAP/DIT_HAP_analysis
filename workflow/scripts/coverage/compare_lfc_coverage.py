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

Bars are drawn by ``cns.barplot`` into the house WIDE panel (150 x 100 pt), so
the Cell palette, typography and bar geometry come from the shared style rather
than from hand-rolled matplotlib calls; the count labels ride on the containers
cns.barplot returns. Every bar carries its absolute LFC-covered gene count. The
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
Version:  1.3.0
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
from typing import Any, cast

# 2. Third-party Imports
import cnsplots as cns
from loguru import logger
from matplotlib.axes import Axes
from matplotlib.container import BarContainer
import numpy as np
import pandas as pd

# 3. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from coverage.core import LFC_TOTAL_LABEL, combine_lfc_coverage_stats
from figures import (
    FIT_ANNOTATION_FONT_PT,
    FURNITURE_COLOR,
    PanelShape,
    apply_house_style,
    fit_panels,
    grid_axes,
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
# Headroom above 100% for the count labels; the plotted data is still the
# untruncated 0-100% range, and the ticks stop at 100.
AXIS_HEADROOM_PCT = 112
AXIS_TICK_VALUES = [0, 20, 40, 60, 80, 100]
X_TICK_ROTATION_DEG = 30
X_AXIS_LABEL = "FYPO viability class"
Y_AXIS_LABEL = "LFC-covered genes (%)"
TITLE = "Gene coverage by FYPO viability class"

# House panel shape: WIDE (150 x 100 pt) for both orientations. The swapped
# chart keeps five class groups legible down the y axis without the empty top
# margin the TALL shape leaves when its aspect is pinned.
PANEL_SHAPE = PanelShape.WIDE

# Legend and footnote positions, in axes fractions below the x axis. The
# vertical chart's rotated tick labels need more room than the swapped one's.
VERTICAL_LEGEND_Y = -0.48
VERTICAL_FOOTNOTE_Y = -0.66
HORIZONTAL_LEGEND_Y = -0.24
HORIZONTAL_FOOTNOTE_Y = -0.42

# The total row's caveat: LD_DIT_HAP's all-gene coverage sits below LD_haploid's,
# so the aggregate bar is starred and the exact counts are footnoted. Both only
# appear when both datasets are plotted.
STAR_DATASET = "LD_DIT_HAP"
REFERENCE_DATASET = "LD_haploid"
STAR_CLASS = LFC_TOTAL_LABEL
CAVEAT_DATASETS = {STAR_DATASET, REFERENCE_DATASET}

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
def _total_lfc_count(comparison: pd.DataFrame, dataset: str) -> int:
    """LFC-covered gene count for one dataset's total row."""
    row = comparison[
        (comparison["dataset"] == dataset) & (comparison["fypo_class"] == LFC_TOTAL_LABEL)
    ]
    return int(row["n_lfc"].iloc[0])


def _plot_frame(
    comparison: pd.DataFrame, datasets: Sequence[str], display: Mapping[str, str]
) -> pd.DataFrame:
    """Long-form frame for cns.barplot: one row per dataset x FYPO class."""
    frame = comparison.loc[
        comparison["dataset"].isin(datasets), ["dataset", "fypo_class", "pct_lfc", "n_lfc"]
    ].copy()
    frame["class_label"] = frame["fypo_class"].map(CLASS_LABELS)
    frame["dataset_label"] = frame["dataset"].map(display)
    frame["coverage_pct"] = frame["pct_lfc"] * PERCENT_PER_UNIT
    return frame


def _count_label(dataset: str, fypo_class: str, count: float, *, show_caveat: bool) -> str:
    """Bar label: the absolute LFC-covered count, starred on the caveat bar."""
    starred = show_caveat and dataset == STAR_DATASET and fypo_class == STAR_CLASS
    return f"{int(count)}{'*' if starred else ''}"


def _draw_count_labels(
    ax: Axes,
    frame: pd.DataFrame,
    classes: list[str],
    datasets: Sequence[str],
    *,
    horizontal: bool,
    show_caveat: bool,
) -> None:
    """Label every cns.barplot container with its absolute LFC-covered count."""
    # seaborn creates one container per hue level, in the hue_order we passed.
    for dataset, container in zip(datasets, ax.containers, strict=True):
        counts = frame[frame["dataset"] == dataset].set_index("fypo_class").loc[classes, "n_lfc"]
        labels = [
            _count_label(dataset, fypo_class, count, show_caveat=show_caveat)
            for fypo_class, count in zip(classes, counts.to_numpy(), strict=True)
        ]
        ax.bar_label(
            # seaborn's categorical plotter only ever appends BarContainers here;
            # the cast narrows Axes.containers' generic Container type.
            cast(BarContainer, container),
            labels=labels,
            padding=1,
            rotation=0 if horizontal else 90,
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
    horizontal: bool = False,
    datasets: Sequence[str] | None = None,
    dataset_labels: Mapping[str, str] | None = None,
) -> None:
    """Draw one grouped LFC-coverage figure in the requested orientation."""
    available = list(dict.fromkeys(comparison["dataset"]))
    selected = list(datasets) if datasets is not None else available
    classes = list(dict.fromkeys(comparison["fypo_class"]))
    display = {dataset: (dataset_labels or {}).get(dataset, dataset) for dataset in selected}
    class_labels = [CLASS_LABELS[fypo_class] for fypo_class in classes]
    colors = series_colors(len(selected))
    accent = house_colors([0])[0]
    show_caveat = CAVEAT_DATASETS.issubset(selected)

    frame = _plot_frame(comparison, selected, display)
    ax = grid_axes(1, 1, labels=[], shape=PANEL_SHAPE)[0]

    bar_kwargs: dict[str, Any] = {
        "data": frame,
        "order": class_labels,
        "hue": "dataset_label",
        "hue_order": [display[dataset] for dataset in selected],
        "palette": colors,
        "ax": ax,
    }

    if horizontal:
        cns.barplot(x="coverage_pct", y="class_label", orient="h", **bar_kwargs)
        _focus_class(ax, classes, accent, horizontal=True)
        ax.set_ylabel(X_AXIS_LABEL)
        ax.set_xlabel(Y_AXIS_LABEL)
        ax.set_xlim(0, AXIS_HEADROOM_PCT)
        ax.set_xticks(AXIS_TICK_VALUES)
        legend_y, footnote_y = HORIZONTAL_LEGEND_Y, HORIZONTAL_FOOTNOTE_Y
    else:
        cns.barplot(x="class_label", y="coverage_pct", orient="v", **bar_kwargs)
        ax.set_xticks(
            np.arange(len(classes)),
            class_labels,
            rotation=X_TICK_ROTATION_DEG,
            ha="right",
        )
        _focus_class(ax, classes, accent, horizontal=False)
        ax.set_xlabel(X_AXIS_LABEL)
        ax.set_ylabel(Y_AXIS_LABEL)
        ax.set_ylim(0, AXIS_HEADROOM_PCT)
        ax.set_yticks(AXIS_TICK_VALUES)
        legend_y, footnote_y = VERTICAL_LEGEND_Y, VERTICAL_FOOTNOTE_Y

    ax.set_title(TITLE)
    _draw_count_labels(
        ax,
        frame,
        classes,
        selected,
        horizontal=horizontal,
        show_caveat=show_caveat,
    )

    handles, legend_labels = ax.get_legend_handles_labels()
    ax.legend(
        handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, legend_y),
        ncol=len(selected),
        frameon=False,
    )

    if show_caveat:
        footnote = (
            f"{STAR_DATASET} {_total_lfc_count(comparison, STAR_DATASET):,} < "
            f"{REFERENCE_DATASET} {_total_lfc_count(comparison, REFERENCE_DATASET):,}"
        )
        ax.text(
            0.5,
            footnote_y,
            footnote,
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=FIT_ANNOTATION_FONT_PT,
            color=FURNITURE_COLOR,
        )

    fit_panels()
    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)


# =============================================================================
# IMPERATIVE SHELL (I/O + orchestration)
# =============================================================================
def read_stats_tables(paths: list[Path]) -> pd.DataFrame:
    """Read and stack the per-dataset stats tables in command-line order."""
    tables = [pd.read_csv(path, sep="\t") for path in paths]
    for path, table in zip(paths, tables, strict=True):
        missing_columns = set(REQUIRED_STATS_COLUMNS) - set(table.columns)
        if missing_columns:
            raise ValueError(f"{path} is missing required columns: {sorted(missing_columns)}")
    return combine_lfc_coverage_stats(tables)


def run(config: CompareLfcCoverageConfig) -> None:
    """Stack the stats tables, write the comparison TSV, then draw the figures."""
    config.validate()
    apply_house_style()

    comparison = read_stats_tables(config.stats)
    comparison.to_csv(config.output_tsv, sep="\t", index=False)
    logger.info(
        f"Wrote {len(comparison)} rows across "
        f"{comparison['dataset'].nunique()} datasets to {config.output_tsv}"
    )

    render_comparison_figure(comparison, config.output_figure.with_suffix(""))
    if config.output_figure_horizontal is not None:
        render_comparison_figure(
            comparison,
            config.output_figure_horizontal.with_suffix(""),
            horizontal=True,
        )
    if config.output_figure_two_libraries is not None:
        render_comparison_figure(
            comparison,
            config.output_figure_two_libraries.with_suffix(""),
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
