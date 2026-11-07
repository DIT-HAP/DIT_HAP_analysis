#!/usr/bin/env python3

"""
STRING Network Enrichment
=========================

The STRING-db half of the non-deterministic network enrichment. Hits the
STRING web API, so its rule stays out of ``rule all`` and every response is
cached under a cache dir — once cached, a re-run needs no network.

Input
-----
- DIT_HAP_all_genes.txt + DIT_HAP_cluster_{n}_genes.txt from the deterministic half

Output
------
- string_enrichment_results.xlsx (one sheet, per-cluster STRING enrichment)

Usage
-----
    python run_string_enrichment.py \\
        --enrichment-dir results/3c_enrichment/HD_DIT_HAP/kmeans_direct9/2026-06-01 \\
        --output-dir .../2026-06-01/network \\
        --cache-dir resources/external/enrichment_cache/HD_DIT_HAP

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-05
Version:  2.0.0
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
# Standard library
import argparse
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

from enrichment.network import (
    StringNetworkConfig,
    run_string_enrichment,
)
from logging_setup import setup_logger

# =============================================================================
# 2. CONSTANTS
# =============================================================================
WORKBOOK_NAME = "string_enrichment_results.xlsx"


# =============================================================================
# 3. CORE LOGIC
# =============================================================================
def write_string_workbook(config: StringNetworkConfig) -> None:
    """Run STRING enrichment and write the results workbook."""
    logger.info("STRING enrichment")
    string_df = run_string_enrichment(config)
    if string_df.empty:
        string_df = pd.DataFrame({"note": ["no STRING results"]})

    output_path = config.output_dir / WORKBOOK_NAME
    with pd.ExcelWriter(output_path) as writer:
        string_df.to_excel(writer, sheet_name="STRING Enrichment", index=False)
    logger.success(f"STRING enrichment complete -> {output_path}")


# =============================================================================
# 4. IMPERATIVE SHELL
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="STRING network enrichment (optional, cached)")
    parser.add_argument(
        "--enrichment-dir", type=Path, required=True, help="Deterministic enrichment results dir"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="Output dir for the STRING workbook"
    )
    parser.add_argument("--cache-dir", type=Path, required=True, help="Cache dir for API responses")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable DEBUG logging")
    return parser.parse_args()


def main() -> int:
    """Validate the boundary once, then run STRING enrichment."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")
    network = StringNetworkConfig(
        enrichment_dir=args.enrichment_dir,
        output_dir=args.output_dir,
        cache_dir=args.cache_dir,
    )

    try:
        network.validate()
        network.output_dir.mkdir(parents=True, exist_ok=True)
        write_string_workbook(network)
    except Exception as error:
        logger.exception(f"STRING enrichment failed: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
