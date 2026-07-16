#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Network Enrichment (STRING + REVIGO)
====================================

The NON-deterministic half of comprehensive_enrichment_analysis.ipynb: STRING-db
functional enrichment and REVIGO semantic-similarity annotation. Both hit external
web APIs, so this rule is optional (not in `rule all`) and caches every response
under a cache dir — once cached, re-runs are deterministic and offline.

Consumes the deterministic enrichment outputs (per-cluster gene lists +
go_enrichment_full.tsv) produced by run_ontology_enrichment.py (Task 4).

Input
-----
- Enrichment output dir from run_ontology_enrichment.py (gene lists + go_enrichment_full.tsv)

Output
------
- string_enrichment_results.xlsx (per-cluster STRING enrichment)
- go_enrichment_full_revigo.tsv (GO enrichment + REVIGO Representative/Eliminated/Dispensability columns)

Usage
-----
    python run_network_enrichment.py \\
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
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.enrichment.pipeline import revigo_analysis, stringdb_enrichment

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
WT_CLUSTER = 9
REVIGO_CUTOFFS = [0.7, 0.5]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class NetworkConfig:
    """Inputs, outputs, and cache location for network enrichment."""
    enrichment_dir: Path
    output_dir: Path
    cache_dir: Path
    wt_cluster: int = WT_CLUSTER

    def validate(self) -> None:
        """Raise ValueError if the upstream enrichment dir is missing."""
        if not self.enrichment_dir.exists():
            raise ValueError(f"Enrichment input dir not found: {self.enrichment_dir}")


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


def _read_gene_list(path: Path) -> list[str]:
    """Read a newline-delimited gene list file into a list, dropping blanks."""
    return [line for line in path.read_text().splitlines() if line.strip()]


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def run_string_enrichment(config: NetworkConfig) -> pd.DataFrame:
    """STRING enrichment for every cluster (background = all genes), concatenated with a Cluster column."""
    bg_genes = _read_gene_list(config.enrichment_dir / "DIT_HAP_all_genes.txt")
    cluster_files = sorted(config.enrichment_dir.glob("DIT_HAP_cluster_*_genes.txt"))

    per_cluster = {}
    for cf in cluster_files:
        cluster = int(cf.stem.replace("DIT_HAP_cluster_", "").replace("_genes", ""))
        genes = _read_gene_list(cf)
        logger.info(f"  STRING cluster {cluster}: {len(genes)} genes")
        per_cluster[cluster] = stringdb_enrichment(genes, bg_genes, cache_dir=config.cache_dir)

    non_empty = {k: v for k, v in per_cluster.items() if not v.empty}
    if not non_empty:
        return pd.DataFrame()
    return pd.concat(non_empty, axis=0).droplevel(1).rename_axis("Cluster").reset_index()


@logger.catch
def annotate_go_with_revigo(config: NetworkConfig) -> pd.DataFrame:
    """Add REVIGO Representative/Eliminated/Dispensability columns to the full GO enrichment table."""
    go_path = config.enrichment_dir / "go_enrichment_full.tsv"
    go_df = pd.read_csv(go_path, sep="\t")
    if go_df.empty:
        return go_df

    annotated_parts = []
    for (cluster, namespace), ns_df in go_df.groupby(["Cluster", "namespace"]):
        ns_df = ns_df.copy()
        for threshold in REVIGO_CUTOFFS:
            revigo = revigo_analysis(ns_df, cut_off=threshold, cache_dir=config.cache_dir)
            revigo = (
                revigo.set_index("Term ID")[["Dispensability", "Eliminated", "Representative"]]
                .add_suffix(f"_{threshold}")
                .reset_index()
            )
            ns_df = ns_df.merge(revigo, left_on="term_id", right_on="Term ID", how="left").drop(columns=["Term ID"])
        annotated_parts.append(ns_df)

    return pd.concat(annotated_parts, axis=0)


@logger.catch
def run_network_enrichment(config: NetworkConfig) -> None:
    """Orchestrate STRING + REVIGO network enrichment, writing results and caching every API response."""
    config.validate()
    config.output_dir.mkdir(parents=True, exist_ok=True)
    config.cache_dir.mkdir(parents=True, exist_ok=True)

    logger.info("STRING enrichment")
    string_df = run_string_enrichment(config)
    with pd.ExcelWriter(config.output_dir / "string_enrichment_results.xlsx") as writer:
        (string_df if not string_df.empty else pd.DataFrame({"note": ["no STRING results"]})).to_excel(
            writer, sheet_name="STRING Enrichment", index=False
        )

    logger.info("REVIGO annotation")
    revigo_df = annotate_go_with_revigo(config)
    revigo_df.to_csv(config.output_dir / "go_enrichment_full_revigo.tsv", sep="\t", index=False)

    logger.success(f"Network enrichment complete -> {config.output_dir}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="STRING + REVIGO network enrichment (optional, cached)")
    parser.add_argument("--enrichment-dir", type=Path, required=True, help="Deterministic enrichment output dir (Task 4)")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output dir for network results")
    parser.add_argument("--cache-dir", type=Path, required=True, help="Cache dir for API responses")
    parser.add_argument("--wt-cluster", type=int, default=WT_CLUSTER, help="WT/background cluster id (default 9)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run network enrichment, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = NetworkConfig(
            enrichment_dir=args.enrichment_dir,
            output_dir=args.output_dir,
            cache_dir=args.cache_dir,
            wt_cluster=args.wt_cluster,
        )
        run_network_enrichment(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
