#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare Coverage Data
======================

Stage 1 of the coverage split: load the insertion-level fitting results +
annotations and the gene-level fitting results, then write two parquet
intermediates consumed by the compute-stats / plot-figures rules:

- annotations.parquet: insertion-level annotations, reindexed onto the
  insertion-level fitting_results' [Chr, Coordinate, Strand, Target] index
  and with duplicate-indexed rows collapsed (see
  workflow.src.coverage.core.resolve_duplicate_annotations) — ready for
  compute_insertion_coverage / compute_per_chromosome_insertion_coverage.
- gene_result.parquet: the full protein-coding gene universe (from the gene
  annotation reference, 1c_annotate.smk) left-joined to the dataset's own gene-level
  fitting results, with legacy um/lam headers normalized to DR/DL.

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
from loguru import logger

# 3. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from io_table import read_parquet, write_parquet  # noqa: E402
from coverage.core import load_insertion_level  # noqa: E402
from release_schema import read_gene_level  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PrepareConfig:
    """Inputs and parquet outputs for the coverage data preparation."""
    fitting_results: Path
    annotations: Path
    gene_level: Path
    annotation_reference: Path
    output_annotations: Path
    output_gene_result: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.fitting_results, self.annotations, self.gene_level, self.annotation_reference]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_annotations, self.output_gene_result]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PrepareConfig) -> None:
    """Load -> build full gene universe from metadata -> left join fitting results -> write parquet intermediates."""
    config.validate()

    gene_result = read_gene_level(config.gene_level)
    _fitting_results, annotations = load_insertion_level(config.fitting_results, config.annotations)

    # gene_level's fitting_results.tsv carries FYPOviability + DeletionLibrary_essentiality,
    # but both are PomBase-era snapshots limited to the genes DIT-HAP happened to cover.
    # Drop them here — FYPOviability and deletion_essentiality come from the annotation
    # reference below (same two sources: PomBase gene metadata + deletion_library_categories
    # .xlsx) applied to the FULL protein-coding gene universe, uncovered genes included.
    # Verified identical to these two columns on the covered subset (2026-07-23
    # coverage-fields verification), and to the reference on the full universe (2026-09-17).
    gene_result = gene_result.drop(columns=["FYPOviability", "DeletionLibrary_essentiality"], errors="ignore")

    # Gene universe + annotation come from the annotation reference (1c_annotate.smk) rather
    # than re-reading PomBase metadata + the deletion-library xlsx here. Only the columns
    # coverage needs are selected, under the reference's own names — so a coverage table
    # reads the same as the reference and the annotated workbook, with no renaming in
    # between. The reference also carries HD_DIT_HAP's own gene-level DR/DL (baked in at
    # build time) and the SGD-derived blocks, none of which belong in a per-dataset
    # coverage analysis.
    reference = read_parquet(config.annotation_reference)
    # protein-only: the reference is built filtered to protein genes, but its experimental
    # blocks are outer-joined, so rows with no current PomBase record carry feature_type=NaN.
    gene_universe = (
        reference.loc[
            reference["feature_type"] == "protein",
            ["gene_name", "product", "characterisation_status", "FYPOviability", "deletion_essentiality"],
        ]
        .rename(columns={"gene_name": "Name"})
        .rename_axis("Systematic ID")
        .reset_index()
    )

    # Left join: all genes from universe, fitting results where available
    gene_result_full = gene_universe.merge(
        gene_result,
        on="Systematic ID",
        how="left",
        suffixes=("_meta", "_fitting")
    )

    # Prefer Name from fitting results if present (it may have been curated), else use metadata
    if "Name_fitting" in gene_result_full.columns:
        gene_result_full["Name"] = gene_result_full["Name_fitting"].fillna(gene_result_full["Name_meta"])
        gene_result_full = gene_result_full.drop(columns=["Name_meta", "Name_fitting"])

    # Genes absent from the deletion library (no deletion-library call was ever made for
    # them) are labeled "Not_determined" rather than left null.
    gene_result_full["deletion_essentiality"] = gene_result_full["deletion_essentiality"].fillna("Not_determined")
    n_not_determined = (gene_result_full["deletion_essentiality"] == "Not_determined").sum()
    logger.info(
        f"Assigned essentiality from the annotation reference: "
        f"{len(gene_result_full) - n_not_determined:,} E/V, {n_not_determined:,} Not_determined"
    )

    logger.info(
        f"Built full gene universe: {len(gene_universe):,} protein-coding genes, "
        f"{gene_result_full['DR'].notna().sum():,} covered (DR not NaN), "
        f"{gene_result_full['DR'].isna().sum():,} not covered"
    )

    if "characterisation_status" in gene_result_full.columns:
        logger.info(f"characterisation_status annotated: {gene_result_full['characterisation_status'].notna().sum():,} genes")
    if "FYPOviability" in gene_result_full.columns:
        n_missing_viability = gene_result_full["FYPOviability"].isna().sum()
        logger.info(f"FYPOviability null count: {n_missing_viability:,} (PomBase's own 'unknown' category covers the rest)")

    write_parquet(annotations, config.output_annotations)
    write_parquet(gene_result_full, config.output_gene_result)

    logger.success(
        f"Prepared coverage data: {len(annotations):,} insertions, {len(gene_result_full):,} genes"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare coverage parquet intermediates")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Insertion-level fitting_results.tsv")
    parser.add_argument("--annotations", type=Path, required=True, help="Insertion-level annotations.tsv(.gz)")
    parser.add_argument("--gene-level", type=Path, required=True, help="Gene-level fitting_results.tsv")
    parser.add_argument("--annotation-reference", type=Path, required=True, help="Gene annotation reference parquet (1c_annotate.smk)")
    parser.add_argument("--output-annotations", type=Path, required=True, help="Output annotations.parquet")
    parser.add_argument("--output-gene-result", type=Path, required=True, help="Output gene_result.parquet")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the preparation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PrepareConfig(
            fitting_results=args.fitting_results,
            annotations=args.annotations,
            gene_level=args.gene_level,
            annotation_reference=args.annotation_reference,
            output_annotations=args.output_annotations,
            output_gene_result=args.output_gene_result,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
