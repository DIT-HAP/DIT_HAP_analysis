#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
REVIGO GO Annotation
=====================

The REVIGO half of the NON-deterministic network enrichment step split off
comprehensive_enrichment_analysis.ipynb (see enrichment_network.smk). Hits the
REVIGO web API, so this rule is optional (not in `rule all`) and caches every
response under a cache dir — once cached, re-runs are deterministic and offline.

Consumes the deterministic enrichment outputs (go_enrichment_full.tsv) produced
by run_ontology_enrichment.py (Task 4).

Input
-----
- Enrichment output dir from run_ontology_enrichment.py (go_enrichment_full.tsv)

Output
------
- go_enrichment_full_revigo.tsv (GO enrichment + REVIGO Representative/Eliminated/Dispensability columns)

Usage
-----
    python run_revigo_annotation.py \\
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

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.enrichment.network import REVIGO_CUTOFFS, WT_CLUSTER, NetworkConfig, annotate_go_with_revigo

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
    """Run REVIGO annotation and write the annotated GO enrichment TSV, caching every API response."""
    config.validate()
    config.output_dir.mkdir(parents=True, exist_ok=True)
    config.cache_dir.mkdir(parents=True, exist_ok=True)

    logger.info("REVIGO annotation")
    revigo_df = annotate_go_with_revigo(config)
    revigo_df.to_csv(config.output_dir / "go_enrichment_full_revigo.tsv", sep="\t", index=False)

    logger.success(f"REVIGO annotation complete -> {config.output_dir}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="REVIGO GO annotation (optional, cached)")
    parser.add_argument("--enrichment-dir", type=Path, required=True, help="Deterministic enrichment output dir (Task 4)")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output dir for REVIGO-annotated results")
    parser.add_argument("--cache-dir", type=Path, required=True, help="Cache dir for API responses")
    parser.add_argument("--wt-cluster", type=int, default=WT_CLUSTER, help="WT/background cluster id (default 9)")
    parser.add_argument(
        "--revigo-cutoffs", type=float, nargs="+", default=list(REVIGO_CUTOFFS),
        help=f"REVIGO semantic-similarity cutoffs, run in order (default {REVIGO_CUTOFFS})",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run REVIGO annotation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = NetworkConfig(
            enrichment_dir=args.enrichment_dir,
            output_dir=args.output_dir,
            cache_dir=args.cache_dir,
            wt_cluster=args.wt_cluster,
            revigo_cutoffs=args.revigo_cutoffs,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
