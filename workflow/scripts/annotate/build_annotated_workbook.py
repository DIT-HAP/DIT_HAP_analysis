#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build Annotated Coverage + Critical Genes Workbook
===================================================

Consolidates coverage detailed_genes.xlsx (16 sheets: All genes + 15 categories)
and verification critical_genes/*.tsv (4 groups) into one annotated Excel workbook.

Output: one master sheet (All genes annotated), 15 category sheets from detailed_genes,
4 critical gene group sheets — 20 sheets total.

Input
-----
- results/2a_coverage/{dataset}/detailed_genes.xlsx
- results/2b_verification/{dataset}/critical_genes/*.tsv
- results/1c_annotation/2026-06-01/2026-08-11/gene_annotation_reference.protein.parquet

Output
------
- {output_xlsx}: consolidated annotated workbook under results/1c_annotation/

Usage
-----
    python build_annotated_workbook.py \\
        --detailed-xlsx results/2a_coverage/HD_DIT_HAP/detailed_genes.xlsx \\
        --critical-dir results/2b_verification/HD_DIT_HAP/critical_genes \\
        --annotation-reference results/1c_annotation/2026-06-01/2026-08-11/gene_annotation_reference.protein.parquet \\
        --output results/1c_annotation/HD_DIT_HAP/HD_DIT_HAP_annotated.xlsx

Author:   Yusheng Yang (guidance) + Claude Opus 5 (implementation)
Date:     2026-08-11
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import argparse
import sys
from pathlib import Path

import pandas as pd
from loguru import logger

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from annotation.core import annotate_table  # noqa: E402
from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# CONFIGURATION
# =============================================================================
EXCEL_SHEET_NAME_MAX_LENGTH = 31


# =============================================================================
# HELPERS
# =============================================================================
def truncate_sheet_name(name: str, used_names: set[str]) -> str:
    """Truncate a sheet name to Excel's 31-character limit, handling collisions."""
    if len(name) <= EXCEL_SHEET_NAME_MAX_LENGTH:
        return name

    base = name[:EXCEL_SHEET_NAME_MAX_LENGTH]
    if base not in used_names:
        return base

    # Handle collision: try base[:28]_01, base[:28]_02, ...
    for i in range(1, 100):
        candidate = f"{name[:28]}_{i:02d}"
        if candidate not in used_names:
            return candidate

    raise ValueError(f"Could not generate a unique truncated name for '{name}' after 99 attempts")


def make_sheet_name_with_count(base_name: str, row_count: int, used_names: set[str]) -> str:
    """Append row count to sheet name, keeping the count inside Excel's 31-char limit.

    The base is trimmed first, so an over-long sheet name loses its tail instead of
    its row count — the count is the part the name exists to convey, and truncating
    from the right turned "essentiality_not_determined (287)" into "(28".
    """
    suffix = f" ({row_count})"
    trimmed_base = base_name[:EXCEL_SHEET_NAME_MAX_LENGTH - len(suffix)]
    return truncate_sheet_name(f"{trimmed_base}{suffix}", used_names)


def annotate_with_reference(
    table: pd.DataFrame,
    reference: pd.DataFrame,
    gene_column: str = "Systematic ID",
) -> pd.DataFrame:
    """Append the reference's columns that `table` does not already carry.

    Anything the table already has came from this same reference — coverage's detailed
    table is built out of it — so re-attaching it would only add identical duplicates,
    and for DR/DL a silently WRONG one: the reference's DR/DL are HD_DIT_HAP's own
    gene-level fits, while a coverage table's belong to the dataset it came from.

    The join itself is annotation.core.annotate_table (one implementation, shared with
    annotate_pombe_genes); genes the reference lacks keep their row with a blank
    annotation rather than raising.

    Places gRNA_DR and gRNA_DL immediately after the input table columns for visibility.
    """
    missing = [c for c in reference.columns if c not in table.columns]
    merged = annotate_table(table, reference, gene_column=gene_column, columns=missing)

    # Move gRNA columns to front of annotation block for visibility
    grna_cols = [c for c in merged.columns if c in ["gRNA_DR", "gRNA_DL"]]
    head = [c for c in merged.columns if c in table.columns]
    tail = [c for c in merged.columns if c not in head and c not in grna_cols]
    return merged[head + grna_cols + tail]


# =============================================================================
# CORE LOGIC
# =============================================================================
def build_workbook(
    detailed_xlsx: Path,
    critical_dir: Path,
    annotation_reference: pd.DataFrame,
    output: Path,
) -> None:
    """Build the consolidated annotated workbook."""
    logger.info(f"Reading detailed genes from {detailed_xlsx}")
    excel_file = pd.ExcelFile(detailed_xlsx)
    sheet_names = excel_file.sheet_names
    logger.info(f"  Found {len(sheet_names)} sheets")

    logger.info(f"Reading critical genes from {critical_dir}")
    critical_files = sorted(critical_dir.glob("*.tsv"))
    logger.info(f"  Found {len(critical_files)} TSV files")

    # Build annotated sheets
    annotated_sheets = {}
    used_names = set()

    # Process detailed_genes sheets (16 sheets)
    for sheet_name in sheet_names:
        df = excel_file.parse(sheet_name)
        annotated = annotate_with_reference(df, annotation_reference)

        output_name = make_sheet_name_with_count(sheet_name, len(df), used_names)
        used_names.add(output_name)
        annotated_sheets[output_name] = annotated
        logger.info(f"  Annotated '{sheet_name}' -> '{output_name}': {len(df):,} rows, {annotated.shape[1]} columns")

    # Process critical_genes files (4 files)
    for tsv_file in critical_files:
        df = pd.read_csv(tsv_file, sep="\t")
        annotated = annotate_with_reference(df, annotation_reference)

        # Use stem as sheet name (e.g., "critical_genes_E2V")
        sheet_name = tsv_file.stem
        output_name = make_sheet_name_with_count(sheet_name, len(df), used_names)
        used_names.add(output_name)
        annotated_sheets[output_name] = annotated
        logger.info(f"  Annotated '{sheet_name}' -> '{output_name}': {len(df):,} rows, {annotated.shape[1]} columns")

    # Write to Excel with "All genes" as the first sheet
    logger.info(f"Writing {len(annotated_sheets)} sheets to {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    # Find the "All genes" sheet (it now has a count suffix)
    all_genes_key = next((k for k in annotated_sheets if k.startswith("All genes")), None)

    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        # Write "All genes" first (master sheet)
        if all_genes_key:
            annotated_sheets[all_genes_key].to_excel(writer, sheet_name=all_genes_key, index=False)
            logger.info(f"  Wrote master sheet '{all_genes_key}'")

        # Write remaining sheets
        for sheet_name, df in annotated_sheets.items():
            if sheet_name != all_genes_key:
                df.to_excel(writer, sheet_name=sheet_name, index=False)

    logger.success(f"Wrote {len(annotated_sheets)} annotated sheets to {output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build annotated workbook from coverage detailed_genes + critical_genes"
    )
    parser.add_argument(
        "--detailed-xlsx",
        type=Path,
        required=True,
        help="Input detailed_genes.xlsx from coverage",
    )
    parser.add_argument(
        "--critical-dir",
        type=Path,
        required=True,
        help="Directory containing critical_genes/*.tsv from verification",
    )
    parser.add_argument(
        "--annotation-reference",
        type=Path,
        required=True,
        help="Annotation reference parquet",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output annotated workbook xlsx",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        # Validate inputs
        for path in [args.detailed_xlsx, args.critical_dir, args.annotation_reference]:
            if not path.exists():
                raise FileNotFoundError(f"Required input does not exist: {path}")

        annotation_reference = read_parquet(args.annotation_reference)
        logger.info(f"Loaded annotation reference: {len(annotation_reference):,} genes × {annotation_reference.shape[1]} columns")

        build_workbook(
            detailed_xlsx=args.detailed_xlsx,
            critical_dir=args.critical_dir,
            annotation_reference=annotation_reference,
            output=args.output,
        )
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
