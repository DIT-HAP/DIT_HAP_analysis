"""
Network Enrichment Core (STRING + REVIGO)
==========================================

Shared config + core logic for the two independent, NON-deterministic network
enrichment steps: STRING-db functional enrichment and REVIGO semantic-similarity
annotation. Both hit external web APIs (via workflow.src.enrichment.pipeline's
stringdb_enrichment / revigo_analysis + cache helpers), so callers should cache
every response under a cache dir — once cached, re-runs are deterministic and
offline.

Consumes the deterministic enrichment outputs (per-cluster gene lists +
go_enrichment_full.tsv) produced by run_ontology_enrichment.py (Task 4).

Usage
-----
    from workflow.src.enrichment.network import (
        NetworkConfig, run_string_enrichment, annotate_go_with_revigo,
    )

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-16
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from dataclasses import dataclass, field
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
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
    revigo_cutoffs: list[float] = field(default_factory=lambda: list(REVIGO_CUTOFFS))

    def validate(self) -> None:
        """Raise ValueError if the upstream enrichment dir is missing."""
        if not self.enrichment_dir.exists():
            raise ValueError(f"Enrichment input dir not found: {self.enrichment_dir}")


# =============================================================================
# HELPERS
# =============================================================================
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
        for threshold in config.revigo_cutoffs:
            revigo = revigo_analysis(ns_df, cut_off=threshold, cache_dir=config.cache_dir)
            revigo = (
                revigo.set_index("Term ID")[["Dispensability", "Eliminated", "Representative"]]
                .add_suffix(f"_{threshold}")
                .reset_index()
            )
            ns_df = ns_df.merge(revigo, left_on="term_id", right_on="Term ID", how="left").drop(columns=["Term ID"])
        annotated_parts.append(ns_df)

    return pd.concat(annotated_parts, axis=0)
