#!/usr/bin/env python3

"""
Cluster Enrichment Dotplot
==========================

The per-variant "face" figure of the enrichment stage: one dot per
cluster x GO term, coloured by -log10(FDR) and sized by the gene ratio. Rows are
the union of each cluster's most significant terms, so the panel answers "which
processes does each cluster stand out for" in one page.

Input
-----
- go_enrichment_full_filtered.tsv (pop_count-capped GO table from finalize_enrichment)

Output
------
- go_enrichment_dotplot.pdf + go_enrichment_dotplot.review.png

Usage
-----
    python plot_enrichment.py \\
        --enrichment .../{variant}/go_enrichment_full_filtered.tsv \\
        --output .../{variant}/go_enrichment_dotplot.pdf \\
        --variant-label kmeans_direct9 \\
        --top-terms-per-cluster 5

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-05
Version:  1.0.0
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
# Standard library
import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

from loguru import logger

# Third-party
import pandas as pd

# Local project (resolved relative to this script, never the current directory)
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR.parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from enrichment.plotting import dotplot_terms, plot_dotplot
from enrichment.schema import ENRICHMENT_COLUMNS
from logging_setup import setup_logger

# =============================================================================
# 2. CONFIGURATION
# =============================================================================


@dataclass(kw_only=True, frozen=True)
class DotplotConfig:
    """Input table, output stem, and the row-selection budget."""

    enrichment: Path
    output: Path
    variant_label: str
    top_terms_per_cluster: int

    @property
    def output_stem(self) -> Path:
        """The extension-free path ``save_dual`` appends .pdf/.review.png to."""
        return self.output.with_suffix("")

    def validate(self) -> None:
        """Raise if the input table is missing, then create the output directory."""
        if not self.enrichment.exists():
            raise ValueError(f"Required input not found: {self.enrichment}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# 3. CORE LOGIC
# =============================================================================
def run(config: DotplotConfig, enrichment: pd.DataFrame) -> None:
    """Select the display terms and render the dotplot."""
    terms = dotplot_terms(enrichment, config.top_terms_per_cluster)
    logger.info(f"{len(terms)} terms across {enrichment[ENRICHMENT_COLUMNS[0]].nunique()} clusters")
    plot_dotplot(enrichment, terms, config.output_stem, title=config.variant_label)


# =============================================================================
# 4. IMPERATIVE SHELL
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Plot one variant's cluster enrichment dotplot")
    parser.add_argument("--enrichment", type=Path, required=True, help="Filtered GO enrichment TSV")
    parser.add_argument("--output", type=Path, required=True, help="Output PDF path")
    parser.add_argument("--variant-label", required=True, help="Variant label shown in the title")
    parser.add_argument(
        "--top-terms-per-cluster",
        type=int,
        required=True,
        help="Terms kept per cluster before the union",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable DEBUG logging")
    return parser.parse_args()


def main() -> int:
    """Validate the boundary once, then render."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")
    config = DotplotConfig(
        enrichment=args.enrichment,
        output=args.output,
        variant_label=args.variant_label,
        top_terms_per_cluster=args.top_terms_per_cluster,
    )

    try:
        config.validate()
        enrichment = pd.read_csv(config.enrichment, sep="\t")
        missing_columns = set(ENRICHMENT_COLUMNS) - set(enrichment.columns)
        if missing_columns:
            logger.error(
                f"Missing required columns in {config.enrichment}: {sorted(missing_columns)}"
            )
            return 1
        run(config, enrichment)
    except Exception as error:
        logger.exception(f"Dotplot failed: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
