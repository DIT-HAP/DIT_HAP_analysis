#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Clustering Distribution — input-space diagnostics
=======================================================

The distribution figure the v2 clustering sweep shipped (design: figure_render/
cluster_dist), rebuilt on the grid-search spine. Six panels ask whether the
space the methods are fit on is one they can plausibly be asked to cluster:

- DR histogram and DL histogram — the one-sided DR tail, the DL point mass at 0;
- the DR/DL cloud coloured by 2D density — the only view where the point mass
  at DL = 0 is visible as a mass rather than an axis tick;
- the standardised coordinates' empirical CDFs on one axis — the only panel
  comparing the two coordinates in SHAPE rather than magnitude;
- one bar per coordinate for its share of total variance — the check on the
  standardisation itself: an axis carrying ~all the variance is ~all the
  distance metric whatever the two are called;
- the same cloud coloured by essentiality call — whether the distinction the
  downstream clusters get judged against is separable here at all.

Raw (DR, DL) is drawn; the standardisation enters through the CDF and the
variance bar.

Input
-----
- annotated_data.parquet + scaled_data.parquet (from prepare_clustering_data)

Output
------
- distribution_stats.tsv: per (axis, scale) descriptive statistics
- distribution.pdf + .review.png (via save_dual)

Usage
-----
    python plot_clustering_distribution.py \\
        --annotated results/clustering/{dataset}/_work/annotated_data.parquet \\
        --scaled results/clustering/{dataset}/_work/scaled_data.parquet \\
        --output-stats results/clustering/{dataset}/distribution_stats.tsv \\
        --output-figure results/clustering/{dataset}/distribution.pdf \\
        --bins 60

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-10-08
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
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import cnsplots as cns  # noqa: E402
from loguru import logger  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from coverage.core import DELETION_ESSENTIALITY_LABELS  # noqa: E402
from figure_render._point_density import point_density  # noqa: E402
from figure_render.scatter import DENSITY_SCATTER_KWS, SCATTERPLOT_KWS  # noqa: E402
from figures import (  # noqa: E402
    FURNITURE_COLOR,
    PanelShape,
    apply_house_style,
    fit_panels,
    grid_axes,
    house_colors,
    panel_labels,
    save_dual,
    series_colors,
)
from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
DR_COLUMN = "DR"
DL_COLUMN = "DL"
DL_DIVISOR = 10  # prepare_clustering_data's DL scaling (dr_clamp is a no-op on the CDF panel: clamped tail)
COORDINATE_COLUMNS = (DR_COLUMN, DL_COLUMN)

# Density markers, enlarged from the shared style: ~4.5k genes pack into a
# 100 px panel, and at the shared s the crowd reads as a wash.
DENSITY_SCATTER_KWS_LARGE = {**DENSITY_SCATTER_KWS, "s": 4}

# The density bar's inset, outside the axes on the right (the cloud's densest
# corner is the DL=0 point mass at high DR — where the default would sit).
DENSITY_CBAR_BOUNDS = (1.06, 0.62, 0.035, 0.30)
DENSITY_CBAR_LABEL = "Density"

# Essentiality panel: raw call -> hue column, draw order, per-meaning colors.
ESSENTIALITY_COLUMN = "RevisedDeletion_essentiality"
ESSENTIALITY_HUE_COLUMN = "Essentiality"
ESSENTIALITY_DRAW_ORDER: tuple[str, ...] = ("V", "Not_determined", "E")

# The stats table's columns, in reading order (v2 distribution lineage).
STAT_COLUMNS: tuple[str, ...] = (
    "axis", "scale", "n", "mean", "std", "min", "q25", "median", "q75", "max", "iqr", "skew", "share_zero",
)
CDF_VALUE_COLUMN = "value"
CDF_CUMULATIVE_COLUMN = "cumulative"


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class DistributionConfig:
    """Spine parquets in; one stats table + one 2x3 diagnostics figure out."""
    annotated: Path
    scaled: Path
    output_stats: Path
    output_figure: Path
    bins: int = 60

    def validate(self) -> None:
        """Raise ValueError if any input is missing, then ensure output dirs exist."""
        for path in [self.annotated, self.scaled]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_stats.parent.mkdir(parents=True, exist_ok=True)
        self.output_figure.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS — statistics
# =============================================================================
def describe_axis(values: pd.Series, *, axis: str, scale: str) -> dict[str, object]:
    """One coordinate's shape descriptors: spread, skew, and the point masses."""
    finite = pd.Series(values).dropna()
    if finite.empty:
        return {"axis": axis, "scale": scale, **{c: np.nan for c in STAT_COLUMNS[2:]}, "n": 0}
    q25, q75 = finite.quantile([0.25, 0.75])
    return {
        "axis": axis,
        "scale": scale,
        "n": int(finite.size),
        "mean": float(finite.mean()),
        "std": float(finite.std()),
        "min": float(finite.min()),
        "q25": float(q25),
        "median": float(finite.median()),
        "q75": float(q75),
        "max": float(finite.max()),
        "iqr": float(q75 - q25),
        "skew": float(finite.skew()),
        "share_zero": float((finite.abs() <= 1e-9).mean()),
    }


def axis_cumulative(values: pd.Series) -> pd.DataFrame:
    """Empirical CDF of one axis: sorted values against their quantile."""
    finite = pd.Series(values).dropna().sort_values()
    return pd.DataFrame(
        {CDF_VALUE_COLUMN: finite.to_numpy(), CDF_CUMULATIVE_COLUMN: np.arange(1, finite.size + 1) / finite.size}
    )


def coordinate_variance_share(frame: pd.DataFrame) -> pd.Series:
    """Each coordinate's share of the total coordinate variance (scaled matrix)."""
    variances = frame[list(COORDINATE_COLUMNS)].var()
    return variances / variances.sum()


# =============================================================================
# HELPERS — panels
# =============================================================================
def render_axis_histogram_panel(ax: Axes, values: pd.Series, *, bins: int, xlabel: str, title: str) -> None:
    """One coordinate's histogram, median marked."""
    finite = pd.Series(values).dropna()
    cns.histplot(data=finite.to_frame("value"), x="value", bins=bins, ax=ax, shrink=0.9)
    median = float(finite.median())
    ax.axvline(median, color=FURNITURE_COLOR, linestyle="--")
    ax.text(
        0.05, 0.95, f"n = {finite.size:,}\nmedian = {median:.3f}\nzero = {(finite.abs() <= 1e-9).mean():.0%}",
        transform=ax.transAxes, va="top",
    )
    ax.set(xlabel=xlabel, ylabel="Genes", title=title)


def render_coordinate_scatter_panel(ax: Axes, frame: pd.DataFrame, *, title: str) -> None:
    """The DR/DL cloud coloured by 2D point density, with a bar in the margin."""
    x = frame[DR_COLUMN].to_numpy(dtype=float)
    y = frame[DL_COLUMN].to_numpy(dtype=float)
    density = point_density(x, y)
    if density is None:
        cns.scatterplot(data=frame, x=DR_COLUMN, y=DL_COLUMN, ax=ax, **SCATTERPLOT_KWS)
    else:
        finite = np.isfinite(density)
        order = np.argsort(density[finite])
        ranked = density[finite][order]
        span = np.ptp(ranked)
        scaled = (ranked - ranked.min()) / span if span > 0 else np.zeros_like(ranked)
        mappable = ax.scatter(
            x[finite][order], y[finite][order], c=scaled, norm=Normalize(vmin=0, vmax=1), **DENSITY_SCATTER_KWS_LARGE
        )
        # scatter.py's helper draws its inset on the cloud's densest corner; for
        # THIS panel the bar belongs in the empty upper-right margin instead
        # (bounds per DENSITY_CBAR_BOUNDS), so draw it here directly.
        cax = ax.inset_axes(DENSITY_CBAR_BOUNDS)
        colorbar = ax.figure.colorbar(mappable, cax=cax, ticks=[0, 1])
        labels = colorbar.ax.set_yticklabels(["low", "high"])
        labels[0].set_verticalalignment("bottom")
        labels[-1].set_verticalalignment("top")
        colorbar.ax.tick_params(length=0, pad=1)
        colorbar.ax.set_title(DENSITY_CBAR_LABEL, loc="center", fontsize=plt.rcParams["ytick.labelsize"])
    ax.set(xlabel=DR_COLUMN, ylabel=DL_COLUMN, title=title)


def render_essentiality_scatter_panel(ax: Axes, frame: pd.DataFrame, *, title: str) -> None:
    """The cloud coloured by essentiality call — is the call separable on these axes at all."""
    labels = frame[ESSENTIALITY_COLUMN].map(lambda call: DELETION_ESSENTIALITY_LABELS.get(call, call))
    hue = frame.assign(**{ESSENTIALITY_HUE_COLUMN: labels})
    observed = set(labels.dropna())
    order = [
        DELETION_ESSENTIALITY_LABELS[value]
        for value in ESSENTIALITY_DRAW_ORDER
        if DELETION_ESSENTIALITY_LABELS[value] in observed
    ]
    order += sorted(str(value) for value in observed - set(order))

    # Meaning-first colours, resolved by palette position rather than order:
    # the essential call takes the warm hit colour (palette index 0), the
    # non-essential bulk the cool baseline (1), missing data the receding grey (7).
    essential = house_colors((0,))[0]
    non_essential = house_colors((1,))[0]
    unknown = house_colors((7,))[0]
    palette = {
        DELETION_ESSENTIALITY_LABELS["E"]: essential,
        DELETION_ESSENTIALITY_LABELS["V"]: non_essential,
        DELETION_ESSENTIALITY_LABELS["Not_determined"]: unknown,
    }
    cns.scatterplot(
        data=hue, x=DR_COLUMN, y=DL_COLUMN, hue=ESSENTIALITY_HUE_COLUMN, hue_order=order,
        palette=palette, ax=ax, **SCATTERPLOT_KWS,
    )
    cns.take_legend_out(ax=ax)
    ax.set(xlabel=DR_COLUMN, ylabel=DL_COLUMN, title=title)


def render_variance_share_panel(ax: Axes, shares: pd.Series, *, title: str) -> None:
    """One bar per coordinate: its share of the total coordinate variance."""
    names = [str(name) for name in shares.index]
    values = shares.to_numpy(dtype=float)
    ax.barh(names, values, color=series_colors(len(names)))
    for index, value in enumerate(values):
        ax.text(value, index, f" {value:.0%}", va="center", ha="left")
    ax.set(xlabel="Share of total variance", ylabel="", title=title)
    ax.set_xlim(0, 1)  # a share on an autoscaled axis overstates the smaller bar
    ax.invert_yaxis()


def render_coordinate_cdf_panel(ax: Axes, frame: pd.DataFrame, series: list[tuple[str, str]]) -> None:
    """Both coordinates' empirical CDFs on one axis, drawn as steps."""
    colors = series_colors(len(series))
    for index, (column, label) in enumerate(series):
        curve = axis_cumulative(frame[column])
        ax.plot(
            curve[CDF_VALUE_COLUMN], curve[CDF_CUMULATIVE_COLUMN],
            drawstyle="steps-post", label=label, color=colors[index],
        )
    ax.legend(loc="best", frameon=False)
    ax.set(xlabel="Standardised coordinate", ylabel="Share of genes at or below", title="Coordinate distributions")


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: DistributionConfig) -> None:
    """Write the per-axis statistics table, then render the 2x3 diagnostics figure."""
    config.validate()

    annotated = read_parquet(config.annotated)
    scaled = read_parquet(config.scaled)
    frame = scaled.join(annotated[[ESSENTIALITY_COLUMN]], how="left")

    stats = pd.DataFrame(
        [
            describe_axis(annotated[DR_COLUMN], axis=DR_COLUMN, scale="raw"),
            describe_axis(annotated[DL_COLUMN], axis=DL_COLUMN, scale="raw"),
            *[
                describe_axis(frame[column], axis=column, scale="scaled")
                for column in COORDINATE_COLUMNS
            ],
        ]
    )[list(STAT_COLUMNS)]
    stats.to_csv(config.output_stats, sep="\t", index=False)

    apply_house_style()
    labels = panel_labels(6)
    axes = grid_axes(2, 3, labels=labels, shape=PanelShape.WIDE)

    render_axis_histogram_panel(axes[0], annotated[DR_COLUMN], bins=config.bins, xlabel=DR_COLUMN, title="Depletion rate")
    render_axis_histogram_panel(axes[1], frame[DL_COLUMN], bins=config.bins, xlabel=DL_COLUMN, title="Depletion lag (scaled)")
    render_coordinate_scatter_panel(axes[2], frame, title="DR vs DL (density)")
    render_coordinate_cdf_panel(
        axes[3], frame,
        [(DR_COLUMN, f"{DR_COLUMN} (clamped)"), (DL_COLUMN, f"{DL_COLUMN} / {DL_DIVISOR}")],
    )
    render_variance_share_panel(axes[4], coordinate_variance_share(scaled), title="Coordinate variance")
    render_essentiality_scatter_panel(axes[5], frame, title="DR vs DL by essentiality")

    fit_panels()  # after all panels drawn: it measures rendered text
    save_dual(config.output_figure.with_suffix(""))
    logger.success(f"distribution -> {config.output_figure} (+ stats {config.output_stats})")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render the clustering input-space diagnostics figure")
    parser.add_argument("--annotated", type=Path, required=True, help="annotated_data.parquet from the spine")
    parser.add_argument("--scaled", type=Path, required=True, help="scaled_data.parquet from the spine")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output distribution_stats.tsv")
    parser.add_argument("--output-figure", type=Path, required=True, help="Output distribution.pdf")
    parser.add_argument("--bins", type=int, default=60, help="Histogram bin count (default 60)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, write stats, render figure."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = DistributionConfig(
            annotated=args.annotated,
            scaled=args.scaled,
            output_stats=args.output_stats,
            output_figure=args.output_figure,
            bins=args.bins,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
