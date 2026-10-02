#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Deletion-Library Paralogues -> Long Table
==========================================

Explodes the Hayles-2013 deletion-library table's `Paralogues` column into one row
per (gene, paralog) pair with both ids resolved to current PomBase systematic ids.

This is not a feature level — nothing merges it into the feature matrix. It is a
rule of its own because the parse is the awkward half (one `|`-joined cell per gene,
the literal `NONE` for "none", a capital-Crick-strand `C` on the 2013 ids) and
because the PAIRS are what a consumer usually wants; `collect_evolutionary_features`
only needs the count, but the coherence attribution's paralog_fraction wants the
gene set, and both would otherwise re-parse the sheet and re-run update_sysIDs().

Input
-----
- The curated deletion-library categories xlsx (Paralogues column)
- PomBase gene_IDs_names_products.tsv, for update_sysIDs()

Output
------
- deletion_library_paralogs.parquet: long table, columns
  gene_systematic_id, gene_name, paralog_systematic_id, paralog_name
  (a name falls back to the systematic id when PomBase has none)

Usage
-----
    python build_deletion_library_paralogs.py \\
        --deletion-library-xlsx resources/curated/deletion_library_categories.xlsx \\
        --gene-meta resources/external/pombase/2026-06-01/Gene_metadata/gene_IDs_names_products.tsv \\
        --output results/1b_features/2026-06-01/deletion_library_paralogs.parquet

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-10-02
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
from io_table import write_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402
from features.assembly import parse_deletion_library_paralogs  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class DeletionLibraryParalogsConfig:
    """Inputs/outputs for the deletion-library paralog long table."""
    deletion_library_xlsx: Path
    gene_meta: Path
    output: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [self.deletion_library_xlsx, self.gene_meta]:
            if not path.exists():
                raise ValueError(f"Required input path does not exist: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: DeletionLibraryParalogsConfig) -> None:
    """Parse the Paralogues column into a long table and write it as parquet."""
    pairs = parse_deletion_library_paralogs(config.deletion_library_xlsx, config.gene_meta)
    write_parquet(pairs, config.output)
    logger.success(
        f"Wrote {len(pairs)} (gene, paralog) pairs over "
        f"{pairs['gene_systematic_id'].nunique()} genes to {config.output}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Explode the deletion-library Paralogues column into a long table")
    parser.add_argument("--deletion-library-xlsx", type=Path, required=True, help="Curated deletion-library categories xlsx")
    parser.add_argument("--gene-meta", type=Path, required=True, help="PomBase gene_IDs_names_products.tsv, for update_sysIDs()")
    parser.add_argument("--output", type=Path, required=True, help="Output long-table parquet")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, write the long table, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = DeletionLibraryParalogsConfig(
            deletion_library_xlsx=args.deletion_library_xlsx,
            gene_meta=args.gene_meta,
            output=args.output,
        )
        config.validate()
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
