#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coherence Annotation Preparer (source -> unified long-table)
============================================================

Stage [1] of the coherence DAG. Reads one PomBase version directory and one
grouping-database source name, dispatches to the matching adapter in
workflow/src/coherence/sources.py, and writes the unified long-table TSV that
every downstream coherence compute/plot stage consumes.

This is a THIN script: all parsing/propagation/reshaping logic lives in the
source adapters (SOURCE_LOADERS). Here we only validate the source name,
dispatch, and write. Adding a database is a one-line change in sources.py +
config.coherence.sources, never here.

Input
-----
- --pombase-dir: a PomBase version directory (the adapter reads its
  ontologies_and_associations/ subtree, e.g. macromolecular_complex_annotation.tsv
  for go_macrocomplex, or go-basic.obo + gene_ontology_annotation.gaf.tsv for
  the go_cc / go_bp GAF-namespace adapters).
- --source: which adapter to run; one of SOURCE_LOADERS keys.

Output
------
- --output: group_annotation_long.tsv with the LONG_TABLE_COLUMNS contract
  (source, group_id, group_name, Systematic ID, Name, n_group_genes) — one row
  per (group, member gene).

Usage
-----
    python prepare_annotation.py \\
        --source go_macrocomplex \\
        --pombase-dir resources/external/pombase/<version> \\
        --output results/coherence/{dataset}/go_macrocomplex/group_annotation_long.tsv

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-23
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from coherence.sources import SOURCE_LOADERS  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# LOGGING SETUP
# =============================================================================
# =============================================================================
# CORE LOGIC
# =============================================================================
def prepare(source: str, pombase_dir: Path) -> pd.DataFrame:
    """Dispatch to the source adapter and return the unified long-table."""
    if source not in SOURCE_LOADERS:
        raise ValueError(f"unknown source {source!r} (have: {sorted(SOURCE_LOADERS)})")
    return SOURCE_LOADERS[source](Path(pombase_dir))


@logger.catch(reraise=True)
def run(source: str, pombase_dir: Path, output: Path) -> None:
    """Prepare the long-table for one source and write it to output."""
    output.parent.mkdir(parents=True, exist_ok=True)
    table = prepare(source, pombase_dir)
    table.to_csv(output, sep="\t", index=False)
    logger.success(f"[{source}] {len(table):,} rows, "
                   f"{table['group_id'].nunique():,} groups -> {output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare a coherence source's unified long-table annotation")
    parser.add_argument("--source", required=True, choices=sorted(SOURCE_LOADERS), help="Grouping-database source adapter to run")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--output", type=Path, required=True, help="Output group_annotation_long.tsv")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: parse args, dispatch to the source adapter, write the long-table."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        run(args.source, args.pombase_dir, args.output)
    except (ValueError, OSError) as e:
        # OSError covers a missing annotation file (the files inside pombase_dir are
        # not individually DAG-tracked) and output write/mkdir failures.
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
