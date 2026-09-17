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

Figures
-------
1. overview              — coverage of every gene, then per essentiality class
2. by_deletion_viability — per class of the curated deletion-library table
3. by_characterisation   — per PomBase characterisation_status
4. insertion_placement   — in-gene vs intergenic, per main chromosome
5. dr / dl_by_essentiality — DR and DL distributions per essentiality class

DR and DL are two files rather than one 3x2 grid because the notebook bins them
on different scales (DR -0.2..1.5 by 0.05, DL 0..15 by 0.5); the grouped
histogram renderer takes one binning per figure.

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-09-17
Version:  4.0.0
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
    DL_BINS,
    DR_BINS,
    IN_GENE_LABEL,
    INTERGENIC_LABEL,
    NOT_COVERED_LABEL,
    characterisation_status_frame,
    coverage_dicts_from_stats_table,
    deletion_viability_frame,
    dr_dl_histogram_frame,
    insertion_placement_frame,
    overall_coverage_frame,
)
from figure_render.composition import render_composition_figure  # noqa: E402
from figure_render.histogram import render_grouped_histogram_figure  # noqa: E402
from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402

# =============================================================================
# CONSTANTS
# =============================================================================
PERCENT_AXIS = "Coverage (%)"


def as_stem(path: Path) -> Path:
    """Drop a trailing .pdf so save_dual appends the suffixes instead of stacking them.

    save_dual takes a stem and writes ``<stem>.pdf`` + ``<stem>.review.png``; handing
    it the rule's declared .pdf output would produce ``<name>.pdf.pdf``.
    """
    return path.with_suffix("") if path.suffix == ".pdf" else path


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PlotFiguresConfig:
    """Inputs (stats TSV + gene_result parquet) and one output stem per figure."""
    stats: Path
    gene_result: Path
    output_overview: Path
    output_deletion_viability: Path
    output_characterisation: Path
    output_insertion_placement: Path
    output_dr_histogram: Path
    output_dl_histogram: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.stats, self.gene_result]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for directory in {path.parent for path in self.output_stems}:
            directory.mkdir(parents=True, exist_ok=True)

    @property
    def output_stems(self) -> list[Path]:
        """Every figure's output stem, in render order."""
        return [
            self.output_overview,
            self.output_deletion_viability,
            self.output_characterisation,
            self.output_insertion_placement,
            self.output_dr_histogram,
            self.output_dl_histogram,
        ]


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotFiguresConfig) -> None:
    """Read stats TSV + gene_result parquet, then render one figure per question."""
    config.validate()

    stats = pd.read_csv(config.stats, sep="\t")
    gene_result = read_parquet(config.gene_result)

    (
        _insertion_coverage,
        gene_coverage,
        _essentiality_coverage,
        per_chromosome,
        characterisation_status_coverage,
        deletion_viability_coverage,
        essentiality_category_coverage,
    ) = coverage_dicts_from_stats_table(stats)

    render_composition_figure(
        overall_coverage_frame(gene_coverage, essentiality_category_coverage),
        config.output_overview,
        category_column="category", percentage_column="covered_pct",
        part_column="covered", whole_column="not_covered",
        part_label=COVERED_LABEL, whole_label=NOT_COVERED_LABEL,
        xlabel=PERCENT_AXIS, ylabel="Essentiality", title="Gene coverage",
    )

    render_composition_figure(
        deletion_viability_frame(deletion_viability_coverage),
        config.output_deletion_viability,
        category_column="category", percentage_column="covered_pct",
        part_column="covered", whole_column="not_covered",
        part_label=COVERED_LABEL, whole_label=NOT_COVERED_LABEL,
        xlabel=PERCENT_AXIS, ylabel="Deletion-library viability",
        title="Gene coverage by deletion viability",
    )

    render_composition_figure(
        characterisation_status_frame(characterisation_status_coverage),
        config.output_characterisation,
        category_column="category", percentage_column="covered_pct",
        part_column="covered", whole_column="not_covered",
        part_label=COVERED_LABEL, whole_label=NOT_COVERED_LABEL,
        xlabel=PERCENT_AXIS, ylabel="Characterisation status",
        title="Gene coverage by characterisation status",
    )

    render_composition_figure(
        insertion_placement_frame(per_chromosome),
        config.output_insertion_placement,
        category_column="category", percentage_column="covered_pct",
        part_column="covered", whole_column="not_covered",
        part_label=IN_GENE_LABEL, whole_label=INTERGENIC_LABEL,
        xlabel="In-gene insertions (%)", ylabel="Chromosome",
        title="Insertion placement by chromosome",
    )

    for feature, bins, stem in (
        ("DR", DR_BINS, config.output_dr_histogram),
        ("DL", DL_BINS, config.output_dl_histogram),
    ):
        render_grouped_histogram_figure(
            dr_dl_histogram_frame(gene_result, feature),
            stem,
            value_column=feature, row_key="stratum", bins=bins,
            xlabel=feature, ylabel="Number of genes",
            # One shared y top: the three strata are meant to be compared, and
            # per-panel scaling makes the 1,144 essential genes' peak look as
            # tall as the 4,513 all-genes peak.
            share_y_range=True,
        )

    logger.success(f"Wrote {len(config.output_stems)} coverage figures")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Plot gene insertion coverage figures")
    parser.add_argument("--stats", type=Path, required=True, help="Input coverage_stats.tsv (composition figures read from here)")
    parser.add_argument("--gene-result", type=Path, required=True, help="Input gene_result.parquet (histograms read from here)")
    parser.add_argument("--output-overview", type=Path, required=True, help="Output stem: overall + per-essentiality coverage")
    parser.add_argument("--output-deletion-viability", type=Path, required=True, help="Output stem: coverage per deletion-library viability")
    parser.add_argument("--output-characterisation", type=Path, required=True, help="Output stem: coverage per characterisation_status")
    parser.add_argument("--output-insertion-placement", type=Path, required=True, help="Output stem: in-gene vs intergenic per chromosome")
    parser.add_argument("--output-dr-histogram", type=Path, required=True, help="Output stem: DR distribution per essentiality")
    parser.add_argument("--output-dl-histogram", type=Path, required=True, help="Output stem: DL distribution per essentiality")
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
            output_overview=as_stem(args.output_overview),
            output_deletion_viability=as_stem(args.output_deletion_viability),
            output_characterisation=as_stem(args.output_characterisation),
            output_insertion_placement=as_stem(args.output_insertion_placement),
            output_dr_histogram=as_stem(args.output_dr_histogram),
            output_dl_histogram=as_stem(args.output_dl_histogram),
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
