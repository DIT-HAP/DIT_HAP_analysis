#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Non-coding RNA Figures
========================

Stage 2b of the non-coding RNA split: read the prepared combined /
nuclear_trnas parquet intermediates and emit the two-page ncRNA figures PDF
(Feature-type donut + tRNA copy-number distribution/DR-by-copy-number
scatter). Depends only on prepare_ncrna_table's output, so it re-runs
independently of the stats rule.

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-22
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

# 2. Third-party Imports
import matplotlib

matplotlib.use("Agg")  # headless: this script only writes a PDF, never displays
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from loguru import logger  # noqa: E402

# 3. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io import read_parquet  # noqa: E402
from workflow.src.noncoding_rna.core import plot_feature_type_donut, plot_trna_summary  # noqa: E402


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class FiguresConfig:
    """Parquet inputs + PDF output for the ncRNA figures."""
    combined: Path
    nuclear_trnas: Path
    output_figures: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [self.combined, self.nuclear_trnas]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_figures.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: FiguresConfig) -> None:
    """Read parquet -> build the 2-page ncRNA figures PDF."""
    config.validate()

    combined = read_parquet(config.combined)
    nuclear_trnas = read_parquet(config.nuclear_trnas)

    fig_donut = plot_feature_type_donut(combined)
    fig_trna = plot_trna_summary(nuclear_trnas)
    with PdfPages(config.output_figures) as pdf:
        pdf.savefig(fig_donut, dpi=300, bbox_inches="tight")
        pdf.savefig(fig_trna, dpi=300, bbox_inches="tight")
    plt.close(fig_donut)
    plt.close(fig_trna)

    logger.success(f"ncRNA figures: {len(combined):,} annotated ncRNA genes, {len(nuclear_trnas):,} nuclear tRNAs")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Non-coding RNA figures (Feature-type donut + tRNA copy-number PDF)")
    parser.add_argument("--combined", type=Path, required=True, help="Input combined.parquet")
    parser.add_argument("--nuclear-trnas", type=Path, required=True, help="Input nuclear_trnas.parquet")
    parser.add_argument("--output-figures", type=Path, required=True, help="Output ncRNA figures PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the figures build, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = FiguresConfig(
            combined=args.combined,
            nuclear_trnas=args.nuclear_trnas,
            output_figures=args.output_figures,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
