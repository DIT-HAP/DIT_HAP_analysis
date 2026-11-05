#!/usr/bin/env python3

"""
Single-Ontology Cluster Enrichment
==================================

Run GO / FYPO / MONDO over-representation (full + slim, goatools FDR-BH) for
every gene cluster, for ONE ontology. The stage fans out per ontology, the
analogue of the clustering stage's per-method fan-out.

Each ontology's DAG and slim association are loaded ONCE and reused across
clusters (the original notebook reloaded them per cluster); enrichment is
deterministic given the same DAG and gene sets, so the results are identical.

Input
-----
- _work/genesets.pkl + _work/id2name.pkl (from prepare_genesets)
- PomBase ontology triples (OBO + GAF/PHAF + slim tables) for the given ontology

Output
------
- {ontology} workbook + nonWT workbook in the results dir
- _work/{ontology}_frames.pkl: (full, slim, nonwt_full, nonwt_slim) frames for finalize

Usage
-----
    python enrich_one_ontology.py \\
        --ontology GO \\
        --genesets results/3c_enrichment/HD_DIT_HAP/kmeans_direct9/2026-06-01/_work/genesets.pkl \\
        --id2name results/3c_enrichment/HD_DIT_HAP/kmeans_direct9/2026-06-01/_work/id2name.pkl \\
        --pombase-dir resources/external/pombase/2026-06-01 \\
        --output-dir results/3c_enrichment/HD_DIT_HAP/kmeans_direct9/2026-06-01 \\
        --work-dir results/3c_enrichment/HD_DIT_HAP/kmeans_direct9/2026-06-01/_work \\
        --wt-cluster 9 --fdr-threshold 0.05

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

from enrichment.cluster_enrichment import (  # noqa: E402
    ONTOLOGIES,
    enrich_all_clusters,
    resolve_ontology,
    write_ontology_workbook,
)
from logging_setup import setup_logger  # noqa: E402

# =============================================================================
# 2. CONFIGURATION
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
    wt_cluster: int
    fdr_threshold: float

    @property
    def ontology_dir(self) -> Path:
        """The PomBase directory holding the OBO/GAF/slim triple."""
        return self.pombase_dir / "ontologies_and_associations"

    def validate(self) -> None:
        """Raise if a required input is missing, then create the output directories."""
        for path in [self.genesets, self.id2name, self.pombase_dir]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.work_dir.mkdir(parents=True, exist_ok=True)


# =============================================================================
# 3. CORE LOGIC
# =============================================================================
def run(config: OntologyConfig) -> None:
    """Enrich every cluster for one ontology, write its workbooks, and pickle the four frames."""
    genesets = pd.read_pickle(config.genesets)
    id2name = pd.read_pickle(config.id2name)

    logger.info(f"{config.ontology} enrichment")
    plan = resolve_ontology(
        config.ontology, config.ontology_dir, config.work_dir, config.fdr_threshold
    )
    full, slim, nonwt_full, nonwt_slim = enrich_all_clusters(
        plan.data,
        genesets.cluster_genes,
        genesets.bg_genes,
        genesets.nonwt_bg_genes,
        plan.load_kwargs,
        plan.enrichment_kwargs,
        {"itemid2name": id2name},
        config.wt_cluster,
    )

    output_dir = config.output_dir
    write_ontology_workbook(
        output_dir / plan.workbook_name, full, slim, plan.full_label, plan.slim_label
    )
    write_ontology_workbook(
        output_dir / plan.nonwt_workbook_name,
        nonwt_full,
        nonwt_slim,
        plan.full_label,
        plan.slim_label,
    )

    pd.to_pickle(
        {"full": full, "slim": slim, "nonwt_full": nonwt_full, "nonwt_slim": nonwt_slim},
        config.work_dir / f"{config.ontology}_frames.pkl",
    )
    logger.success(f"[{config.ontology}] enrichment complete -> {config.output_dir}")


# =============================================================================
# 4. IMPERATIVE SHELL
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run one ontology's cluster enrichment (deterministic)"
    )
    parser.add_argument("--ontology", required=True, choices=ONTOLOGIES, help="Ontology to enrich")
    parser.add_argument("--genesets", type=Path, required=True, help="Pickled ClusterGeneSets")
    parser.add_argument("--id2name", type=Path, required=True, help="Pickled id->name dict")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="Results dir for the workbooks"
    )
    parser.add_argument(
        "--work-dir", type=Path, required=True, help="Work dir for GAF intermediates + frames"
    )
    parser.add_argument("--wt-cluster", type=int, required=True, help="WT/background cluster id")
    parser.add_argument("--fdr-threshold", type=float, required=True, help="FDR-BH alpha")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable DEBUG logging")
    return parser.parse_args()


def main() -> int:
    """Validate the boundary once, then hand the trusted config to the core."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")
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

    try:
        config.validate()
        run(config)
    except Exception as error:
        logger.exception(f"[{config.ontology}] enrichment failed: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
