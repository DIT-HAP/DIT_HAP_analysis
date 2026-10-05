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
- gene_result.parquet: the full protein-coding gene universe with its annotation
  columns (all from the gene annotation reference, 1c_annotate.smk) left-joined to
  the dataset's own gene-level DR/DL, so uncovered genes survive as DR=NaN rows.
  Legacy um/lam headers are normalized to DR/DL on the way in.

Input
-----
- {release_dir}/insertion_level/fitting_results.tsv and annotations.tsv.gz
- {release_dir}/gene_level/fitting_results.tsv — this dataset's own DR/DL
- results/1c_annotation/{pombase_version}/{sgd_version}/gene_annotation_reference.protein.parquet

Output
------
- results/2a_coverage/{dataset}/_work/annotations.parquet
- results/2a_coverage/{dataset}/_work/gene_result.parquet

Usage
-----
    python workflow/scripts/coverage/prepare_coverage_data.py \\
        --fitting-results .../insertion_level/fitting_results.tsv \\
        --annotations .../insertion_level/annotations.tsv.gz \\
        --gene-level .../gene_level/fitting_results.tsv \\
        --annotation-reference results/1c_annotation/2026-06-01/2026-08-11/gene_annotation_reference.protein.parquet \\
        --output-annotations results/2a_coverage/HD_DIT_HAP/_work/annotations.parquet \\
        --output-gene-result results/2a_coverage/HD_DIT_HAP/_work/gene_result.parquet

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
@dataclass(kw_only=True, slots=True, frozen=True)
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
    """Load -> take the gene universe from the annotation reference -> left join this dataset's DR/DL -> write parquet."""
    config.validate()

    # DR/DL come from the dataset's own gene-level table — the annotation reference's
    # are dataset-independent (HD_DIT_HAP's), so taking them here would silently
    # mislabel every other dataset (DR correlates with HD_DIT_HAP at r = 0.61-0.95 on
    # the other three released datasets). Nothing else is read from it: Name matches the
    # reference's gene_name on all four released datasets, and its FYPOviability /
    # DeletionLibrary_essentiality are a covered-genes-only snapshot superseded by the
    # reference's full-universe versions.
    gene_result = read_gene_level(config.gene_level)[["Systematic ID", "DR", "DL"]]
    _fitting_results, annotations = load_insertion_level(config.fitting_results, config.annotations)

    # Gene universe + annotation come from the annotation reference (1c_annotate.smk) — every
    # column it carries, under its own names, so a coverage table reads the same as the
    # reference and the annotated workbook, with no renaming in between. The reference's
    # HD_DIT_HAP_DR/DL stay (named for the dataset they came from, so they sit beside this
    # dataset's own DR/DL as a comparison, not a collision).
    reference = read_parquet(config.annotation_reference)
    # protein-only: the reference is built filtered to protein genes, but its experimental
    # blocks are outer-joined, so rows with no current PomBase record carry feature_type=NaN.
    gene_result_full = (
        reference.loc[reference["feature_type"] == "protein"]
        .rename(columns={"gene_name": "Name"})
        .rename_axis("Systematic ID")
        .reset_index()
        .merge(gene_result, on="Systematic ID", how="left")
    )

    # Genes absent from the deletion library (no deletion-library call was ever made for
    # them) are labeled "Not_determined" rather than left null.
    gene_result_full["deletion_essentiality"] = gene_result_full["deletion_essentiality"].fillna("Not_determined")
    n_not_determined = (gene_result_full["deletion_essentiality"] == "Not_determined").sum()
    logger.info(
        f"Assigned essentiality from the annotation reference: "
        f"{len(gene_result_full) - n_not_determined:,} E/V, {n_not_determined:,} Not_determined"
    )

    logger.info(
        f"Built full gene universe: {len(gene_result_full):,} protein-coding genes, "
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
