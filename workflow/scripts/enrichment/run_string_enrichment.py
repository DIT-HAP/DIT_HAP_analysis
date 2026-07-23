#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
STRING Network Enrichment
==========================

The STRING-db half of the NON-deterministic network enrichment step split off
comprehensive_enrichment_analysis.ipynb (see enrichment_network.smk). Hits the
STRING-db web API, so this rule is optional (not in `rule all`) and caches every
response under a cache dir — once cached, re-runs are deterministic and offline.

Consumes the deterministic enrichment outputs (per-cluster gene lists +
DIT_HAP_all_genes.txt) produced by run_ontology_enrichment.py (Task 4).

Input
-----
- Enrichment output dir from run_ontology_enrichment.py (per-cluster gene lists + all-genes list)

Output
------
- string_enrichment_results.xlsx (per-cluster STRING enrichment)

Usage
-----
    python run_string_enrichment.py \\
        --enrichment-dir results/enrichment/raw/{dataset}/{pombase_version} \\
        --output-dir results/enrichment/network/{dataset}/{pombase_version} \\
        --cache-dir resources/external/enrichment_cache/{dataset}

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-16
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
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.enrichment.network import WT_CLUSTER, NetworkConfig, run_string_enrichment

# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def run(config: NetworkConfig) -> None:
    """Run STRING enrichment and write the results workbook, caching every API response."""
    config.validate()
    config.output_dir.mkdir(parents=True, exist_ok=True)
    config.cache_dir.mkdir(parents=True, exist_ok=True)

    logger.info("STRING enrichment")
    string_df = run_string_enrichment(config)
    with pd.ExcelWriter(config.output_dir / "string_enrichment_results.xlsx") as writer:
        (string_df if not string_df.empty else pd.DataFrame({"note": ["no STRING results"]})).to_excel(
            writer, sheet_name="STRING Enrichment", index=False
        )

    logger.success(f"STRING enrichment complete -> {config.output_dir}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="STRING network enrichment (optional, cached)")
    parser.add_argument("--enrichment-dir", type=Path, required=True, help="Deterministic enrichment output dir (Task 4)")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output dir for STRING results")
    parser.add_argument("--cache-dir", type=Path, required=True, help="Cache dir for API responses")
    parser.add_argument("--wt-cluster", type=int, default=WT_CLUSTER, help="WT/background cluster id (default 9)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run STRING enrichment, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = NetworkConfig(
            enrichment_dir=args.enrichment_dir,
            output_dir=args.output_dir,
            cache_dir=args.cache_dir,
            wt_cluster=args.wt_cluster,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
