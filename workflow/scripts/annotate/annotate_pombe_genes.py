#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pombe Gene Table Annotation
=============================

Appends annotation columns to any table that carries a column of pombe
systematic IDs, to help judge what a gene does and whether it is essential:
budding yeast ortholog name and null-mutant essentiality, pombe FYPO and
deletion-library essentiality, GO-slim terms, complex membership, PFAM domains.

Gene IDs are matched verbatim against the annotation reference — no synonym
resolution is attempted, so the input is expected to carry current systematic
IDs. Unmatched IDs are reported rather than silently dropped, since a wrong
--gene-column or a non-coding gene slipping in looks exactly like a table with
no annotations.

Input
-----
- Any tsv/csv/xlsx table with a column of pombe systematic IDs
- gene_annotation_reference.parquet (built by build_annotation_reference.py)

Output
------
- The input table with annotation columns appended (tsv/csv/xlsx by extension)

Usage
-----
    python annotate_pombe_genes.py \\
        --input my_gene_list.tsv \\
        --gene-column gene_systematic_id \\
        --output my_gene_list.annotated.tsv

    # only some annotation columns, and drop rows that have no annotation
    python annotate_pombe_genes.py \\
        --input results/clustering/HD_DIT_HAP/direct/final_clusters.tsv \\
        --gene-column Gene_id \\
        --columns gene_name Sc_ortholog_name Sc_essentiality \\
        --drop-unmatched \\
        --output annotated_clusters.tsv

Author:   Yusheng Yang (guidance) + Claude Opus 5 (implementation)
Date:     2026-08-11
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
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from annotation.core import annotate_table, summarise_match  # noqa: E402
from io_table import read_file, read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# How many unmatched gene ids to name in the log before truncating.
_MAX_REPORTED_UNMATCHED = 20


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class AnnotateConfig:
    """Inputs/outputs and join options for table annotation."""

    input_table: Path
    gene_column: str
    annotation_reference: Path
    output: Path
    columns: list[str] | None
    drop_unmatched: bool

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [self.input_table, self.annotation_reference]:
            if not path.exists():
                raise ValueError(f"Required input path does not exist: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def write_table(table: pd.DataFrame, output: Path) -> None:
    """Write a table as tsv/csv/xlsx, dispatching on the output extension."""
    match output.suffix.lower():
        case ".tsv" | ".txt":
            table.to_csv(output, sep="\t", index=False)
        case ".csv":
            table.to_csv(output, index=False)
        case ".xlsx":
            table.to_excel(output, index=False)
        case _:
            raise ValueError(f"Unsupported output extension: {output.suffix} (use .tsv/.csv/.xlsx)")


def report_match(table: pd.DataFrame, reference: pd.DataFrame, gene_column: str) -> None:
    """Log how many rows matched the reference and which gene ids did not."""
    matched, unmatched = summarise_match(table, reference, gene_column=gene_column)
    logger.info(f"Matched {matched:,}/{len(table):,} rows against the annotation reference")

    duplicates = int(table[gene_column].duplicated().sum())
    if duplicates:
        logger.info(f"{duplicates:,} rows repeat a gene id; each is annotated in place")

    if not unmatched:
        return

    shown = ", ".join(unmatched[:_MAX_REPORTED_UNMATCHED])
    suffix = f" (+{len(unmatched) - _MAX_REPORTED_UNMATCHED} more)" if len(unmatched) > _MAX_REPORTED_UNMATCHED else ""
    logger.warning(f"{len(unmatched):,} gene id(s) had no annotation: {shown}{suffix}")
    logger.warning("Check --gene-column, and whether these are current systematic IDs of coding genes")


# =============================================================================
# CORE LOGIC
# =============================================================================
def run(config: AnnotateConfig) -> None:
    """Join annotation columns onto the input table and write the annotated result.

    Deliberately not wrapped in @logger.catch: the expected failures here (wrong
    --gene-column, absent annotation column) already carry actionable messages, and
    main() reports them. A full traceback would bury the message that matters.
    """
    table = read_file(config.input_table)
    reference = read_parquet(config.annotation_reference)
    logger.info(
        f"Read {len(table):,} rows from {config.input_table} and "
        f"{len(reference):,} annotated genes from {config.annotation_reference}"
    )

    report_match(table, reference, config.gene_column)

    annotated = annotate_table(
        table,
        reference,
        gene_column=config.gene_column,
        columns=config.columns,
        drop_unmatched=config.drop_unmatched,
    )
    added = annotated.shape[1] - table.shape[1]
    write_table(annotated, config.output)
    logger.success(f"Wrote {len(annotated):,} rows x {added} added columns to {config.output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Append annotation columns to a table of pombe genes")
    parser.add_argument("--input", type=Path, required=True, dest="input_table", help="Input table (tsv/csv/xlsx)")
    parser.add_argument(
        "--gene-column", required=True, help="Name of the column holding pombe systematic IDs"
    )
    parser.add_argument(
        "--annotation-reference",
        type=Path,
        required=True,
        help="Annotation reference parquet (see build_annotation_reference.py)",
    )
    parser.add_argument("--output", type=Path, required=True, help="Output annotated table (tsv/csv/xlsx)")
    parser.add_argument(
        "--columns", nargs="+", default=None, help="Only append these annotation columns (default: all)"
    )
    parser.add_argument(
        "--drop-unmatched",
        action="store_true",
        help="Drop rows with no annotation (default: keep them with empty annotation columns)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, annotate the table, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = AnnotateConfig(
            input_table=args.input_table,
            gene_column=args.gene_column,
            annotation_reference=args.annotation_reference,
            output=args.output,
            columns=args.columns,
            drop_unmatched=args.drop_unmatched,
        )
        config.validate()
        run(config)
    except (ValueError, KeyError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
