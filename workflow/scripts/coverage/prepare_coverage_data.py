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
- gene_result.parquet: gene-level fitting results, with legacy um/lam headers
  normalized to DR/DL.

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
import pandas as pd

# 3. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from io_table import read_parquet, write_parquet  # noqa: E402
from coverage.core import load_gene_level, load_insertion_level  # noqa: E402
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
    gene_metadata: Path
    deletion_library_xlsx: Path
    output_annotations: Path
    output_gene_result: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.fitting_results, self.annotations, self.gene_level, self.gene_metadata, self.deletion_library_xlsx]:
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

    gene_result = load_gene_level(config.gene_level)
    _fitting_results, annotations = load_insertion_level(config.fitting_results, config.annotations)

    # gene_level's fitting_results.tsv carries FYPOviability + DeletionLibrary_essentiality,
    # but both are PomBase-era snapshots limited to the genes DIT-HAP happened to cover.
    # Drop them here — deletion_viability (from gene_metadata, below) and essentiality
    # (from deletion_library_categories.xlsx, below) are the same two facts sourced
    # directly and applied to the FULL protein-coding gene universe, uncovered genes
    # included. Verified byte-for-byte equal to these two on the covered subset (see
    # 2026-07-23 coverage-fields verification): FYPOviability == deletion_viability after
    # normalizing the "condition-dependent"/"depends_on_conditions" label spelling; native
    # DeletionLibrary_essentiality == deletion_library_categories.xlsx's "Gene
    # dispensability. This study", with every native "Not_determined" exactly matching a
    # gene absent from that xlsx.
    gene_result = gene_result.drop(columns=["FYPOviability", "DeletionLibrary_essentiality"], errors="ignore")

    # Load full protein-coding gene universe from gene metadata
    gene_metadata = read_parquet(config.gene_metadata)
    protein_genes = gene_metadata[gene_metadata["feature_type"] == "protein"].copy()

    # Start with full protein-coding gene list, then left join fitting results
    # This ensures uncovered genes (no DR/DL) are present as DR=NaN rows
    full_gene_cols = ["systematic_id", "name", "characterisation_status", "deletion_viability"]
    available_cols = [c for c in full_gene_cols if c in protein_genes.columns]
    gene_universe = protein_genes[available_cols].copy()
    gene_universe = gene_universe.rename(columns={"systematic_id": "Systematic ID", "name": "Name"})

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

    # essentiality: sourced straight from deletion_library_categories.xlsx (the curated
    # Hayles-derived dispensability study) for every gene in the universe, not just the
    # ones DIT-HAP's curve fitting covered. Genes absent from the xlsx (no deletion-library
    # call was ever made for them) are labeled "Not_determined" rather than left null.
    deletion_library = pd.read_excel(config.deletion_library_xlsx)
    dl_essentiality = deletion_library.set_index("Systematic ID")["Gene dispensability. This study"]
    gene_result_full["essentiality"] = (
        gene_result_full["Systematic ID"].map(dl_essentiality).fillna("Not_determined")
    )
    n_not_determined = (gene_result_full["essentiality"] == "Not_determined").sum()
    logger.info(
        f"Assigned essentiality from deletion_library_categories.xlsx: "
        f"{len(gene_result_full) - n_not_determined:,} E/V, {n_not_determined:,} Not_determined"
    )

    logger.info(
        f"Built full gene universe: {len(gene_universe):,} protein-coding genes, "
        f"{gene_result_full['DR'].notna().sum():,} covered (DR not NaN), "
        f"{gene_result_full['DR'].isna().sum():,} not covered"
    )

    if "characterisation_status" in gene_result_full.columns:
        logger.info(f"characterisation_status annotated: {gene_result_full['characterisation_status'].notna().sum():,} genes")
    if "deletion_viability" in gene_result_full.columns:
        n_missing_viability = gene_result_full["deletion_viability"].isna().sum()
        logger.info(f"deletion_viability null count: {n_missing_viability:,} (PomBase's own 'unknown' category covers the rest)")

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
    parser.add_argument("--gene-metadata", type=Path, required=True, help="Gene metadata parquet with characterisation_status")
    parser.add_argument("--deletion-library-xlsx", type=Path, required=True, help="Deletion library categories Excel file")
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
            gene_metadata=args.gene_metadata,
            deletion_library_xlsx=args.deletion_library_xlsx,
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
