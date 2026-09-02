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
- results/coverage/{dataset}/detailed_genes.xlsx
- results/verification/{dataset}/critical_genes/*.tsv
- gene_annotation_reference.parquet

Output
------
- {output_xlsx}: consolidated annotated workbook under results/annotation/

Usage
-----
    python build_annotated_workbook.py \\
        --detailed-xlsx results/coverage/HD_DIT_HAP/detailed_genes.xlsx \\
        --critical-dir results/verification/HD_DIT_HAP/critical_genes \\
        --annotation-reference results/annotation/2026-06-01/2026-08-11/gene_annotation_reference.parquet \\
        --output results/annotation/HD_DIT_HAP/HD_DIT_HAP_annotated.xlsx

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

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io_table import read_parquet


# =============================================================================
# CONFIGURATION
# =============================================================================
EXCEL_SHEET_NAME_MAX_LENGTH = 31


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )


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
    """Append row count to sheet name and truncate if needed to fit Excel's 31-char limit."""
    name_with_count = f"{base_name} ({row_count})"
    return truncate_sheet_name(name_with_count, used_names)


def annotate_with_suffix_on_collision(
    table: pd.DataFrame,
    reference: pd.DataFrame,
    gene_column: str = "Systematic ID",
) -> pd.DataFrame:
    """Annotate a table, adding _annotation suffix to annotation columns that collide with existing ones.

    Places gRNA_DR and gRNA_DL immediately after the input table columns for visibility.
    """
    annotation = reference.loc[table[gene_column], :]
    annotation = annotation.reset_index(drop=True)

    # Detect column collisions
    collisions = set(table.columns) & set(annotation.columns)
    if collisions:
        rename_map = {col: f"{col}_annotation" for col in collisions}
        annotation = annotation.rename(columns=rename_map)

    # Move gRNA columns to front of annotation block for visibility
    grna_cols = [c for c in annotation.columns if c in ["gRNA_DR", "gRNA_DL"]]
    other_cols = [c for c in annotation.columns if c not in grna_cols]
    annotation = annotation[grna_cols + other_cols]

    return pd.concat([table.reset_index(drop=True), annotation], axis=1)


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
        annotated = annotate_with_suffix_on_collision(df, annotation_reference)

        output_name = make_sheet_name_with_count(sheet_name, len(df), used_names)
        used_names.add(output_name)
        annotated_sheets[output_name] = annotated
        logger.info(f"  Annotated '{sheet_name}' -> '{output_name}': {len(df):,} rows, {annotated.shape[1]} columns")

    # Process critical_genes files (4 files)
    for tsv_file in critical_files:
        df = pd.read_csv(tsv_file, sep="\t")
        annotated = annotate_with_suffix_on_collision(df, annotation_reference)

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
