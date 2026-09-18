#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Coverage Figures
======================

Stage 2b of the coverage split: render the coverage figures with cnsplots, one
PDF per question. Every breakdown is a composition figure (a percentage bar plus
one part/whole donut per category) read straight from coverage_stats.tsv, so a
figure can never disagree with the numbers compute_coverage_stats wrote. The
DR/DL histograms are the exception: they need the per-gene values that the
aggregated stats table does not carry, so they read the gene_result parquet.

Each column in coverage.core.DIMENSION_LABELS gets the SAME two figures, broken
down the same way, so no annotation dimension is analysed differently from the
others:

  coverage_by_{column}.pdf      — composition: "All genes", then one bar per value
  coverage_dr_by_{column}.pdf   — DR distribution, one panel per value
  coverage_dl_by_{column}.pdf   — DL distribution, one panel per value

plus one figure that is not per-gene-dimension (coverage_insertion_placement.pdf,
in-gene vs intergenic per main chromosome). File names are derived from the column
names, so the figure set follows DIMENSION_LABELS instead of being listed here.

DR and DL are two files rather than one grid because the notebook bins them on
different scales (DR -1.45..0.2 by 0.05, DL 0..15 by 0.5); the grouped histogram
renderer takes one binning per figure.

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-09-17
Version:  5.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Third-party Imports
import pandas as pd
from loguru import logger

# 3. Local Imports (relative path resolution)
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from coverage.core import (  # noqa: E402
    COVERED_LABEL,
    DIMENSION_LABELS,
    DL_BINS,
    DR_BINS,
    IN_GENE_LABEL,
    INTERGENIC_LABEL,
    NOT_COVERED_LABEL,
    coverage_dicts_from_stats_table,
    dimension_coverage_frame,
    dr_dl_histogram_frame,
    insertion_placement_frame,
)
from figure_render.composition import render_composition_figure  # noqa: E402
from figure_render.histogram import render_grouped_histogram_figure  # noqa: E402
from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402

# =============================================================================
# CONSTANTS
# =============================================================================
PERCENT_AXIS = "Coverage (%)"

# Display (ylabel, title) for each dimension's composition figure, keyed by the
# gene_result column. A dimension missing here falls back to its column name
# rather than failing, so adding one to DIMENSION_LABELS still renders.
DIMENSION_AXES = {
    "characterisation_status": ("Characterisation status", "Gene coverage by characterisation status"),
    "FYPOviability": ("FYPO viability", "Gene coverage by FYPO viability"),
    "deletion_essentiality": ("Deletion essentiality", "Gene coverage by deletion essentiality"),
}

HISTOGRAM_FEATURES = (("DR", DR_BINS), ("DL", DL_BINS))


def figure_stems(output_dir: Path) -> dict[str, Path]:
    """Every figure's output stem: a composition + two DR/DL histograms per dimension, plus insertion placement."""
    stems = {
        f"composition_{column}": output_dir / f"coverage_by_{column}"
        for column in DIMENSION_LABELS
    }
    stems |= {
        f"histogram_{feature.lower()}_{column}": output_dir / f"coverage_{feature.lower()}_by_{column}"
        for column in DIMENSION_LABELS
        for feature, _ in HISTOGRAM_FEATURES
    }
    stems["insertion_placement"] = output_dir / "coverage_insertion_placement"
    return stems


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PlotFiguresConfig:
    """Inputs (stats TSV + gene_result parquet) and the directory the figures are written to."""
    stats: Path
    gene_result: Path
    output_dir: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [self.stats, self.gene_result]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_dir.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotFiguresConfig) -> None:
    """Read stats TSV + gene_result parquet, then render one figure per question."""
    config.validate()

    stats = pd.read_csv(config.stats, sep="\t")
    gene_result = read_parquet(config.gene_result)
    stems = figure_stems(config.output_dir)

    gene_coverage, per_chromosome, dimension_coverage = coverage_dicts_from_stats_table(stats)

    for column in DIMENSION_LABELS:
        ylabel, title = DIMENSION_AXES.get(column, (column, f"Gene coverage by {column}"))

        render_composition_figure(
            dimension_coverage_frame(gene_coverage, dimension_coverage[column], DIMENSION_LABELS[column]),
            stems[f"composition_{column}"],
            category_column="category", percentage_column="covered_pct",
            part_column="covered", whole_column="not_covered",
            part_label=COVERED_LABEL, whole_label=NOT_COVERED_LABEL,
            xlabel=PERCENT_AXIS, ylabel=ylabel, title=title,
        )

        for feature, bins in HISTOGRAM_FEATURES:
            render_grouped_histogram_figure(
                dr_dl_histogram_frame(gene_result, feature, column, DIMENSION_LABELS[column]),
                stems[f"histogram_{feature.lower()}_{column}"],
                value_column=feature, row_key="stratum", bins=bins,
                xlabel=feature, ylabel="Number of genes",
                # Per-panel y, NOT shared: dimension values differ by an order of
                # magnitude (characterisation_status runs 2,456 down to 11), and one
                # shared top flattens everything but the biggest panel to nothing.
                # The x axis stays shared, so the shapes are still comparable.
                share_y_range=False,
            )

    render_composition_figure(
        insertion_placement_frame(per_chromosome),
        stems["insertion_placement"],
        category_column="category", percentage_column="covered_pct",
        part_column="covered", whole_column="not_covered",
        part_label=IN_GENE_LABEL, whole_label=INTERGENIC_LABEL,
        xlabel="In-gene insertions (%)", ylabel="Chromosome",
        title="Insertion placement by chromosome",
    )

    logger.success(f"Wrote {len(stems)} coverage figures")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Plot gene insertion coverage figures")
    parser.add_argument("--stats", type=Path, required=True, help="Input coverage_stats.tsv (composition figures read from here)")
    parser.add_argument("--gene-result", type=Path, required=True, help="Input gene_result.parquet (histograms read from here)")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory the figures are written to")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render every figure, report the outcome."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotFiguresConfig(
            stats=args.stats,
            gene_result=args.gene_result,
            output_dir=args.output_dir,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
