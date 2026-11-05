#!/usr/bin/env python3

"""
Enrichment Preprocessing (the spine)
====================================

First step of the GO/FYPO/MONDO enrichment stage: turn the curated
``final_clusters.tsv`` into background / nonWT / per-cluster gene sets, build
the id->name map, and write the gene-list text files. Emits the two pickles the
per-ontology jobs read, so the ontology fan-out needs no clustering knowledge.

``final_clusters.tsv`` is sourced through ``final_clusters_path(dataset,
variant)`` in the clustering stage; every finalize variant builds one.

Input
-----
- final_clusters.tsv (Systematic ID, cluster, ...) from a finalize variant
- PomBase gene metadata + the curated deletion-library workbook

Output
------
- DIT_HAP_all_genes.txt, DIT_HAP_cluster_{n}_genes.txt, DIT_HAP_cluster_genes_matrix.txt
- _work/genesets.pkl (ClusterGeneSets), _work/id2name.pkl (dict)

Usage
-----
    python prepare_genesets.py \\
        --final-clusters results/3b_clustering/HD_DIT_HAP/kmeans_direct9/final_clusters.tsv \\
        --pombase-dir resources/external/pombase/2026-06-01 \\
        --deletion-library-xlsx resources/curated/deletion_library_categories.xlsx \\
        --output-dir results/3c_enrichment/HD_DIT_HAP/kmeans_direct9/2026-06-01 \\
        --work-dir results/3c_enrichment/HD_DIT_HAP/kmeans_direct9/2026-06-01/_work \\
        --wt-cluster 9

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-05
Version:  2.0.0
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
# Standard library
import argparse
from dataclasses import dataclass
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

from enrichment.cluster_enrichment import (
    CLUSTER_COLUMN,
    load_cluster_genesets,
    write_gene_lists,
)
from enrichment.ontology import GeneMetaConfig
from enrichment.schema import FINAL_CLUSTER_COLUMNS
from io_table import read_file
from logging_setup import setup_logger

# =============================================================================
# 2. CONFIGURATION
# =============================================================================


@dataclass(kw_only=True, frozen=True)
class PrepareConfig:
    """Inputs, output/work dirs, and gene-set parameters for enrichment preprocessing."""

    final_clusters: Path
    pombase_dir: Path
    deletion_library_xlsx: Path
    output_dir: Path
    work_dir: Path
    wt_cluster: int
    cluster_column: str = CLUSTER_COLUMN

    def validate(self) -> None:
        """Raise if a required input is missing, then create the output directories."""
        for path in [self.final_clusters, self.pombase_dir, self.deletion_library_xlsx]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.work_dir.mkdir(parents=True, exist_ok=True)


# =============================================================================
# 3. CORE LOGIC
# =============================================================================
def run(config: PrepareConfig, clusters: pd.DataFrame) -> None:
    """Build gene sets + id2name, write the gene lists, and pickle both for the ontology jobs."""
    gene_meta = GeneMetaConfig(
        gene_ids_names_products=config.pombase_dir
        / "Gene_metadata"
        / "gene_IDs_names_products.tsv",
        deletion_library_essentiality=config.deletion_library_xlsx,
    ).load_data()

    genesets = load_cluster_genesets(clusters, config.cluster_column, config.wt_cluster)
    write_gene_lists(genesets.cluster_genes, genesets.bg_genes, config.output_dir)
    logger.info(
        f"{len(genesets.cluster_genes)} clusters, {len(genesets.bg_genes)} background genes"
    )

    pd.to_pickle(genesets, config.work_dir / "genesets.pkl")
    pd.to_pickle(gene_meta.id2name, config.work_dir / "id2name.pkl")
    logger.success(f"Prepared gene sets + id2name -> {config.work_dir}")


# =============================================================================
# 4. IMPERATIVE SHELL
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Prepare enrichment gene sets (spine)")
    parser.add_argument(
        "--final-clusters", type=Path, required=True, help="Curated final_clusters.tsv"
    )
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument(
        "--deletion-library-xlsx",
        type=Path,
        required=True,
        help="Curated deletion library workbook",
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="Results dir for the gene lists"
    )
    parser.add_argument("--work-dir", type=Path, required=True, help="Work dir for the pickles")
    parser.add_argument("--wt-cluster", type=int, required=True, help="WT/background cluster id")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable DEBUG logging")
    return parser.parse_args()


def main() -> int:
    """Validate the boundary once, then hand the trusted frame to the core."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")
    config = PrepareConfig(
        final_clusters=args.final_clusters,
        pombase_dir=args.pombase_dir,
        deletion_library_xlsx=args.deletion_library_xlsx,
        output_dir=args.output_dir,
        work_dir=args.work_dir,
        wt_cluster=args.wt_cluster,
    )

    try:
        config.validate()
        clusters = read_file(config.final_clusters)
        missing_columns = set(FINAL_CLUSTER_COLUMNS) - set(clusters.columns)
        if missing_columns:
            logger.error(
                f"Missing required columns in {config.final_clusters}: {sorted(missing_columns)}"
            )
            return 1
        run(config, clusters)
    except Exception as error:
        logger.exception(f"Enrichment preprocessing failed: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
