#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coherence Interactive Scatter — an Altair HTML explorer
=======================================================

The static group_scatter.pdf answers "where do this term's genes sit in DR-DL
space?" for a namelist fixed at config time. Answering it for a term you did not
think of means editing config and re-running the pipeline, and the static figure
cannot say anything about a gene you point at.

This is the same picture as an interactive document: the genome-wide DR/DL cloud
is drawn once, a dropdown picks the term, and the term's genes are highlighted on
top of it. Hovering any point — background or highlighted — shows that gene's
fitness record. The output is a single self-contained HTML file.

The term list is every group in the source's metrics table, most coherent first,
so the dropdown opens on the strongest signal. Genes are the term's *scored*
members (the DR<threshold set the coherence z-score was computed on), joined to
the upstream fitting table for the per-gene detail shown on hover.

Input
-----
- coherence_metrics.parquet: one source's per-group metrics (source, group_id,
  group_name, n_scored_members, median_pairwise_distance_z, q_value).
- group_annotation_long.tsv: that source's prepared group -> member long table.
- fitting_results.tsv: the upstream per-gene fitness table (see coherence/io.py).

Output
------
- interactive_scatter.html: one self-contained page. It loads vega-embed from a
  CDN, so opening it needs a network connection the first time.

Usage
-----
    python plot_interactive_scatter.py \\
        --metrics results/3a_coherence/{dataset}/{source}/coherence_metrics.parquet \\
        --annotation results/3a_coherence/{dataset}/{source}/group_annotation_long.tsv \\
        --fitting-results .../gene_level/fitting_results.tsv \\
        --source go_macrocomplex \\
        --output tmp/{source}/interactive_scatter.html

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-09-22
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
import altair as alt
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from coherence.io import load_fitting_results, load_long_table  # noqa: E402
from figures import house_colors  # noqa: E402
from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# The genome cloud is context, not data: drawn at the same marker size as the
# highlighted genes — a smaller one made the cloud read as a different kind of
# object rather than as the same measurement in context — but in a light grey and
# translucent enough that a few thousand overlapping markers stay a background.
_BACKGROUND_LIGHT_GREY = "#c9c9c9"
_BACKGROUND_OPACITY = 0.35

# The highlighted genes are the subject of the figure, and a term is 3-300 genes,
# so they can afford to be large. The cloud matches this size.
_MEMBER_SIZE = 70
_MEMBER_OPACITY = 0.85
_BACKGROUND_SIZE = _MEMBER_SIZE

# Floor for the DL axis, matching the attribution figure's. Most genes sit at
# exactly DL = 0, so an axis ending at 0 draws half of those markers outside the
# frame; DR has no such pile-up and keeps its own limits.
_DL_FLOOR = -0.05

# Chart footprint in pixels. Larger than a journal panel: this is read on screen
# and every point has to be big enough to hover.
_CHART_WIDTH = 620
_CHART_HEIGHT = 460

_REQUIRED_METRIC_COLUMNS = ["group_name", "median_pairwise_distance_z", "q_value"]

# Per-gene columns carried into the tooltip, in display order. Every one comes
# from the upstream fitting table except the two normalized coordinates.
_GENE_TOOLTIP = [
    alt.Tooltip("Name:N", title="Gene"),
    alt.Tooltip("Systematic ID:N", title="Systematic ID"),
    alt.Tooltip("DR:Q", title="DR", format=".3f"),
    alt.Tooltip("DL:Q", title="DL", format=".3f"),
    alt.Tooltip("norm_DR:Q", title="norm DR", format=".3f"),
    alt.Tooltip("norm_DL:Q", title="norm DL/10", format=".3f"),
    alt.Tooltip("R2:Q", title="R²", format=".3f"),
    alt.Tooltip("FYPOviability:N", title="Viability"),
]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Inputs and output for the interactive scatter page."""
    metrics: Path
    annotation: Path
    fitting_results: Path
    source: str
    output: Path

    def validate(self) -> None:
        """Raise ValueError if an input is missing, then create the output dir."""
        for path in [self.metrics, self.annotation, self.fitting_results]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
def member_table(
    metrics: pd.DataFrame, long_table: pd.DataFrame, fitting: pd.DataFrame, source: str
) -> pd.DataFrame:
    """The term -> gene rows to draw: one row per (group, scored member)."""
    # Inner join on `Systematic ID`, so a member the fitting table could not fit
    # is dropped rather than drawn at a fabricated coordinate — the same members
    # the coherence z-score was computed on.
    members = long_table[long_table["source"] == source].copy()
    detail = fitting[[
        "Systematic ID", "DR", "DL", "norm_DR", "norm_DL", "R2", "FYPOviability",
    ]]
    joined = members.merge(detail, on="Systematic ID", how="inner")
    return joined.merge(
        metrics[["group_id", "median_pairwise_distance_z", "q_value", "n_scored_members"]],
        on="group_id", how="inner",
    )


def term_options(metrics: pd.DataFrame) -> list[str]:
    """Dropdown entries: every group, most coherent first."""
    ordered = metrics.sort_values("median_pairwise_distance_z")
    return [str(name) for name in ordered["group_name"].dropna().unique()]


def build_chart(background: pd.DataFrame, members: pd.DataFrame, options: list[str]) -> alt.Chart:
    """The layered page: genome cloud under the term-picked genes."""
    # The term picker is a selection on `group_name`, used as a TRANSFORM FILTER on
    # the highlighted layer. Filtering the layer rather than pre-slicing the frame
    # is what keeps all the terms in one page, so the dropdown a reader changes
    # needs no re-render.
    picker = alt.selection_point(
        fields=["group_name"],
        bind=alt.binding_select(options=options, name="Term  "),
        value=options[0],
    )

    cloud = (
        alt.Chart(background)
        .mark_circle(size=_BACKGROUND_SIZE, color=_BACKGROUND_LIGHT_GREY, opacity=_BACKGROUND_OPACITY)
        .encode(
            x=alt.X("norm_DR:Q", title="norm DR"),
            y=alt.Y("norm_DL:Q", title="norm DL/10", scale=alt.Scale(domainMin=_DL_FLOOR)),
            tooltip=_GENE_TOOLTIP,
        )
    )
    highlighted = (
        alt.Chart(members)
        .transform_filter(picker)
        .mark_circle(size=_MEMBER_SIZE, color=house_colors((0,))[0], opacity=_MEMBER_OPACITY)
        .encode(
            x=alt.X("norm_DR:Q", title="norm DR"),
            y=alt.Y("norm_DL:Q", title="norm DL/10", scale=alt.Scale(domainMin=_DL_FLOOR)),
            tooltip=_GENE_TOOLTIP + [
                alt.Tooltip("group_name:N", title="Term"),
                alt.Tooltip("median_pairwise_distance_z:Q", title="z-score", format=".2f"),
                alt.Tooltip("q_value:Q", title="q (BH)", format=".3g"),
            ],
        )
    )
    return (
        (cloud + highlighted)
        .add_params(picker)
        .properties(
            width=_CHART_WIDTH,
            height=_CHART_HEIGHT,
            title="Coherence members in DR–DL space (hover a point for its gene record)",
        )
        .interactive()
    )


@logger.catch(reraise=True)
def run(config: PlotConfig) -> None:
    """Load -> build the terms and gene rows -> save the interactive page."""
    config.validate()
    metrics = read_parquet(config.metrics)
    missing = [column for column in _REQUIRED_METRIC_COLUMNS if column not in metrics.columns]
    if missing:
        raise ValueError(f"metrics table missing required column(s) {missing} (have: {list(metrics.columns)})")

    long_table = load_long_table(config.annotation)
    fitting = load_fitting_results(config.fitting_results)

    members = member_table(metrics, long_table, fitting, config.source)
    options = term_options(metrics)
    if not options:
        raise ValueError(f"no groups for source '{config.source}' in {config.metrics}")

    logger.info(
        f"{len(metrics):,} groups in '{config.source}', {len(members):,} member rows, "
        f"{len(fitting):,} background genes; opening on '{options[0]}'"
    )
    chart = build_chart(fitting, members, options)
    chart.save(config.output)
    logger.success(f"Wrote {config.output} ({len(members):,} member rows over {len(options):,} terms)")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Interactive Altair scatter of a source's coherence terms")
    parser.add_argument("--metrics", type=Path, required=True, help="A source's coherence_metrics.parquet")
    parser.add_argument("--annotation", type=Path, required=True, help="That source's group_annotation_long.tsv")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Upstream fitting_results.tsv")
    parser.add_argument("--source", type=str, required=True, help="Source name (e.g. go_macrocomplex)")
    parser.add_argument("--output", type=Path, required=True, help="Output HTML page")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render the page, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotConfig(
            metrics=args.metrics,
            annotation=args.annotation,
            fitting_results=args.fitting_results,
            source=args.source,
            output=args.output,
        )
        run(config)
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
