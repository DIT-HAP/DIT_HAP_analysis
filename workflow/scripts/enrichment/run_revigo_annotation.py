#!/usr/bin/env python3

"""
REVIGO GO Annotation
====================

The REVIGO half of the non-deterministic network enrichment. Hits the REVIGO
web API, so its rule stays out of ``rule all`` and every response is cached
under a cache dir — once cached, a re-run needs no network.

Input
-----
- go_enrichment_full.tsv from the deterministic half

Output
------
- go_enrichment_full_revigo.tsv (GO enrichment + one REVIGO
  Representative/Eliminated/Dispensability triplet per configured cutoff)

Usage
-----
    python run_revigo_annotation.py \\
        --enrichment-dir results/3c_enrichment/HD_DIT_HAP/kmeans_direct9/2026-06-01 \\
        --output-dir .../2026-06-01/network \\
        --cache-dir resources/external/enrichment_cache/HD_DIT_HAP \\
        --revigo-cutoffs 0.7 0.5

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

# Third-party
from loguru import logger

# Local project (resolved relative to this script, never the current directory)
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR.parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from enrichment.network import (
    RevigoNetworkConfig,
    annotate_go_with_revigo,
)
from logging_setup import setup_logger

# =============================================================================
# 2. CONSTANTS
# =============================================================================
OUTPUT_NAME = "go_enrichment_full_revigo.tsv"


# =============================================================================
# 3. CORE LOGIC
# =============================================================================
def write_revigo_table(config: RevigoNetworkConfig) -> None:
    """Run REVIGO annotation and write the annotated GO table."""
    logger.info("REVIGO annotation")
    revigo_df = annotate_go_with_revigo(config)
    output_path = config.output_dir / OUTPUT_NAME
    revigo_df.to_csv(output_path, sep="\t", index=False)
    logger.success(f"REVIGO annotation complete -> {output_path}")


# =============================================================================
# 4. IMPERATIVE SHELL
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="REVIGO GO annotation (optional, cached)")
    parser.add_argument(
        "--enrichment-dir", type=Path, required=True, help="Deterministic enrichment results dir"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="Output dir for the REVIGO table"
    )
    parser.add_argument("--cache-dir", type=Path, required=True, help="Cache dir for API responses")
    parser.add_argument(
        "--revigo-cutoffs",
        type=float,
        nargs="+",
        required=True,
        help="REVIGO semantic-similarity cutoffs, run in order",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable DEBUG logging")
    return parser.parse_args()


def main() -> int:
    """Validate the boundary once, then run REVIGO annotation."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")
    network = RevigoNetworkConfig(
        enrichment_dir=args.enrichment_dir,
        output_dir=args.output_dir,
        cache_dir=args.cache_dir,
        revigo_cutoffs=args.revigo_cutoffs,
    )

    try:
        network.validate()
        network.output_dir.mkdir(parents=True, exist_ok=True)
        write_revigo_table(network)
    except Exception as error:
        logger.exception(f"REVIGO annotation failed: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
