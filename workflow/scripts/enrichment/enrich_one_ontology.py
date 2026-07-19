#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Single-Ontology Cluster Enrichment
=====================================

Runs GO / FYPO / MONDO over-representation analysis (full + slim, goatools
FDR-BH) for every gene cluster, for ONE ontology. Fanned out per ontology by
enrichment.smk (design doc §5) — the analogue of clustering.smk's method split.

Each ontology's DAG + slim association is loaded ONCE and reused across clusters
(the source notebook reloaded per cluster); enrichment is deterministic given
the same DAG and gene sets, so results are identical.

Input
-----
- _work/genesets.pkl + _work/id2name.pkl (from prepare_genesets)
- PomBase ontology triples (OBO + GAF/PHAF + slim tables) for the given ontology

Output
------
- {ontology} workbook + nonWT workbook (canonical names) in the raw dir
- _work/{ontology}_frames.pkl: (full, slim, nonwt_full, nonwt_slim) frames for finalize

Usage
-----
    python enrich_one_ontology.py \\
        --ontology GO \\
        --genesets .../_work/genesets.pkl --id2name .../_work/id2name.pkl \\
        --pombase-dir resources/external/pombase/2025-10-01 \\
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
from workflow.src.enrichment.cluster_enrichment import (
    FDR_THRESHOLD,
    ONTOLOGIES,
    WT_CLUSTER,
    enrich_all_clusters,
    resolve_ontology,
    write_ontology_workbook,
)


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class OntologyConfig:
    """Inputs, output/work dirs, and parameters for a single ontology's enrichment."""
    ontology: str
    genesets: Path
    id2name: Path
    pombase_dir: Path
    output_dir: Path
    work_dir: Path
    wt_cluster: int = WT_CLUSTER
    fdr_threshold: float = FDR_THRESHOLD

    @property
    def ontology_dir(self) -> Path:
        return self.pombase_dir / "ontologies_and_associations"

    def validate(self) -> None:
        """Raise ValueError on an unknown ontology or missing input; ensure output/work dirs exist."""
        if self.ontology not in ONTOLOGIES:
            raise ValueError(f"Unknown ontology: {self.ontology!r} (expected one of {ONTOLOGIES})")
        for path in [self.genesets, self.id2name, self.pombase_dir]:
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
def run(config: OntologyConfig) -> None:
    """Enrich all clusters for one ontology, write its workbooks, and pickle the 4 frames for finalize."""
    config.validate()
    genesets = pd.read_pickle(config.genesets)
    id2name = pd.read_pickle(config.id2name)
    format_kwargs = {"itemid2name": id2name}

    logger.info(f"{config.ontology} enrichment")
    plan = resolve_ontology(config.ontology, config.ontology_dir, config.work_dir, config.fdr_threshold)
    full, slim, nonwt_full, nonwt_slim = enrich_all_clusters(
        plan.data, genesets.cluster_genes, genesets.bg_genes, genesets.nonwt_bg_genes,
        plan.load_kwargs, plan.enrichment_kwargs, format_kwargs, config.wt_cluster,
    )

    write_ontology_workbook(config.output_dir / plan.workbook_name, full, slim, plan.full_label, plan.slim_label)
    write_ontology_workbook(config.output_dir / plan.nonwt_workbook_name, nonwt_full, nonwt_slim, plan.full_label, plan.slim_label)

    pd.to_pickle(
        {"full": full, "slim": slim, "nonwt_full": nonwt_full, "nonwt_slim": nonwt_slim},
        config.work_dir / f"{config.ontology}_frames.pkl",
    )
    logger.success(f"[{config.ontology}] enrichment complete -> {config.output_dir}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Run one ontology's cluster enrichment (deterministic)")
    parser.add_argument("--ontology", required=True, choices=ONTOLOGIES, help="Ontology to enrich")
    parser.add_argument("--genesets", type=Path, required=True, help="Pickled ClusterGeneSets (from prepare)")
    parser.add_argument("--id2name", type=Path, required=True, help="Pickled id->name dict (from prepare)")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--output-dir", type=Path, required=True, help="Raw results dir (workbooks)")
    parser.add_argument("--work-dir", type=Path, required=True, help="Work dir for GAF intermediates + frame pickle")
    parser.add_argument("--wt-cluster", type=int, default=WT_CLUSTER, help="WT/background cluster id (default 9)")
    parser.add_argument("--fdr-threshold", type=float, default=FDR_THRESHOLD, help="FDR-BH alpha (default 0.05)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run one ontology's enrichment, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = OntologyConfig(
            ontology=args.ontology,
            genesets=args.genesets,
            id2name=args.id2name,
            pombase_dir=args.pombase_dir,
            output_dir=args.output_dir,
            work_dir=args.work_dir,
            wt_cluster=args.wt_cluster,
            fdr_threshold=args.fdr_threshold,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
