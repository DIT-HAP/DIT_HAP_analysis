#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enrichment Preprocessing (the spine)
======================================

First stage of the split GO/FYPO/MONDO enrichment pipeline: reads the curated
final_clusters.tsv into background / nonWT / per-cluster gene sets, builds the
id->name map, and writes the gene-list txt files. Emits pickles for the three
per-ontology enrichment jobs (design doc §5).

final_clusters.tsv is an UN-BUILDABLE human-curated input (design doc §8): if
missing, Snakemake reports "missing input" — run the finalize notebook first.

Input
-----
- Curated final_clusters.tsv (Systematic ID, cluster, ...)
- PomBase gene metadata + curated deletion-library table (for id->name mapping)

Output
------
- gene-list txt files (DIT_HAP_all_genes.txt, per-cluster, matrix) in the raw dir
- _work/genesets.pkl (ClusterGeneSets), _work/id2name.pkl (dict)

Usage
-----
    python prepare_genesets.py \\
        --final-clusters resources/curated/final_clusters.tsv \\
        --pombase-dir resources/external/pombase/2025-10-01 \\
        --deletion-library-xlsx resources/curated/deletion_library_categories.xlsx \\
        --output-dir results/enrichment/raw/{dataset}/{version} \\
        --work-dir results/enrichment/raw/{dataset}/{version}/_work

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-17
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
from workflow.src.enrichment.cluster_enrichment import CLUSTER_COLUMN, WT_CLUSTER, load_cluster_genesets, write_gene_lists
from workflow.src.enrichment.ontology import GeneMetaConfig


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PrepareConfig:
    """Inputs, output/work dirs, and gene-set parameters for enrichment preprocessing."""
    final_clusters: Path
    pombase_dir: Path
    deletion_library_xlsx: Path
    output_dir: Path
    work_dir: Path
    cluster_column: str = CLUSTER_COLUMN
    wt_cluster: int = WT_CLUSTER

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output/work dirs exist."""
        for path in [self.final_clusters, self.pombase_dir, self.deletion_library_xlsx]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.work_dir.mkdir(parents=True, exist_ok=True)


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
@logger.catch(reraise=True)
def run(config: PrepareConfig) -> None:
    """Load gene sets + id->name map, write gene lists, and pickle both for the ontology jobs."""
    config.validate()

    gene_meta = GeneMetaConfig(
        gene_IDs_names_products=config.pombase_dir / "Gene_metadata" / "gene_IDs_names_products.tsv",
        deletion_library_essentiality=config.deletion_library_xlsx,
    ).load_data()

    genesets = load_cluster_genesets(config.final_clusters, config.cluster_column, config.wt_cluster)
    write_gene_lists(genesets.cluster_genes, genesets.bg_genes, config.output_dir)
    logger.info(f"{len(genesets.cluster_genes)} clusters, {len(genesets.bg_genes)} background genes")

    pd.to_pickle(genesets, config.work_dir / "genesets.pkl")
    pd.to_pickle(gene_meta.id2name, config.work_dir / "id2name.pkl")
    logger.success(f"Prepared gene sets + id2name -> {config.work_dir}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare enrichment gene sets (spine)")
    parser.add_argument("--final-clusters", type=Path, required=True, help="Curated final_clusters.tsv")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--deletion-library-xlsx", type=Path, required=True, help="Curated deletion library categories xlsx")
    parser.add_argument("--output-dir", type=Path, required=True, help="Raw results dir (gene lists)")
    parser.add_argument("--work-dir", type=Path, required=True, help="Work dir for pickled genesets/id2name")
    parser.add_argument("--wt-cluster", type=int, default=WT_CLUSTER, help="WT/background cluster id (default 9)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run preprocessing, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PrepareConfig(
            final_clusters=args.final_clusters,
            pombase_dir=args.pombase_dir,
            deletion_library_xlsx=args.deletion_library_xlsx,
            output_dir=args.output_dir,
            work_dir=args.work_dir,
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
