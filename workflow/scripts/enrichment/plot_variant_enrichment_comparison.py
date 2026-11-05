#!/usr/bin/env python3

"""
Cross-Variant Enrichment Comparison
===================================

Two single-page figures that answer "does the clustering variant change the
enrichment story?":

- a term x variant heatmap of the strongest -log10(FDR) each variant reaches for
  a term (rows = union of each variant's top terms; blank = not significant)
- grouped bars of the significant-term count per cluster for every variant

Input
-----
- one go_enrichment_full_filtered.tsv per variant, plus a label per file

Output
------
- variant_enrichment_heatmap.pdf + .review.png
- variant_enrichment_counts.pdf + .review.png

Usage
-----
    python plot_variant_enrichment_comparison.py \\
        --input .../{variant}/go_enrichment_full_filtered.tsv \\
                .../{variant}/go_enrichment_full_filtered.tsv \\
        --label kmeans_direct9 hc_direct9 \\
        --output-dir results/3c_enrichment/HD_DIT_HAP \\
        --top-terms 15 --fdr-threshold 0.05

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

from enrichment.plotting import (  # noqa: E402
    plot_variant_counts,
    plot_variant_heatmap,
    significant_term_counts,
    term_display_names,
    variant_heatmap_matrix,
)
from enrichment.schema import ENRICHMENT_COLUMNS  # noqa: E402
from logging_setup import setup_logger  # noqa: E402

# =============================================================================
# 2. CONSTANTS
# =============================================================================
HEATMAP_STEM = "variant_enrichment_heatmap"
COUNTS_STEM = "variant_enrichment_counts"


# =============================================================================
# 3. CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class ComparisonConfig:
    """Per-variant input tables, their labels, and the comparison parameters."""

    inputs: list[Path]
    labels: list[str]
    output_dir: Path
    top_terms: int
    fdr_threshold: float

    def validate(self) -> None:
        """Raise on a label/input length mismatch or a missing input, then create the output dir."""
        if len(self.inputs) != len(self.labels):
            raise ValueError(f"Got {len(self.inputs)} inputs but {len(self.labels)} labels")
        for path in self.inputs:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_dir.mkdir(parents=True, exist_ok=True)


# =============================================================================
# 4. CORE LOGIC
# =============================================================================
def run(config: ComparisonConfig, per_variant: dict[str, pd.DataFrame]) -> None:
    """Render the heatmap and the significant-term counts."""
    matrix = variant_heatmap_matrix(per_variant, config.top_terms)
    plot_variant_heatmap(matrix, term_display_names(per_variant), config.output_dir / HEATMAP_STEM)

    counts = significant_term_counts(per_variant, config.fdr_threshold)
    plot_variant_counts(counts, config.output_dir / COUNTS_STEM, config.fdr_threshold)


# =============================================================================
# 5. IMPERATIVE SHELL
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Compare GO enrichment across clustering variants")
    parser.add_argument(
        "--input", type=Path, nargs="+", required=True, help="Filtered GO TSVs, one per variant"
    )
    parser.add_argument(
        "--label", nargs="+", required=True, help="Variant label per input, same order"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="Directory for the two figures"
    )
    parser.add_argument(
        "--top-terms", type=int, required=True, help="Terms kept per variant before the union"
    )
    parser.add_argument(
        "--fdr-threshold", type=float, required=True, help="FDR cut for the count panel"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable DEBUG logging")
    return parser.parse_args()


def _load_variants(config: ComparisonConfig) -> dict[str, pd.DataFrame]:
    """Read every variant table and validate its columns at the boundary."""
    per_variant: dict[str, pd.DataFrame] = {}
    for path, label in zip(config.inputs, config.labels, strict=True):
        frame = pd.read_csv(path, sep="\t")
        missing_columns = set(ENRICHMENT_COLUMNS) - set(frame.columns)
        if missing_columns:
            raise ValueError(f"Missing required columns in {path}: {sorted(missing_columns)}")
        per_variant[label] = frame
    return per_variant


def main() -> int:
    """Validate the boundary once, then render both comparison figures."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")
    config = ComparisonConfig(
        inputs=args.input,
        labels=args.label,
        output_dir=args.output_dir,
        top_terms=args.top_terms,
        fdr_threshold=args.fdr_threshold,
    )

    try:
        config.validate()
        run(config, _load_variants(config))
    except Exception as error:
        logger.exception(f"Variant comparison failed: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
