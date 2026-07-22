#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Network-Level Feature Collection
==================================

Assembles GO term richness (from the GO GAF) and BioGrid PPI/GI degree per
coding gene. Reads the coding-gene set from the DNA-level parquet and loads the
GO DAG/GAF to derive gene2go.

Input
-----
- A PomBase version directory (go-basic.obo, gene_ontology_annotation.gaf.tsv, slim tables)
- A BioGrid interaction table
- DNA-level features parquet (for the coding-gene set)

Output
------
- network_features.parquet: per-gene network feature table (indexed by gene id)

Usage
-----
    python collect_network_features.py \\
        --pombase-dir resources/external/pombase/2025-10-01 \\
        --biogrid-tsv resources/external/biogrid/BIOGRID-....tab3.txt \\
        --dna-features results/features/2025-10-01/_levels/dna_features.parquet \\
        --output results/features/2025-10-01/_levels/network_features.parquet

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

# 2. Third-party Imports
from loguru import logger

# 3. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.enrichment.ontology import OntologyDataConfig, load_ontology_data
from workflow.src.features.assembly import collect_network_level_features, read_coding_genes


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class NetworkConfig:
    """Inputs/outputs for network-level feature collection."""
    pombase_dir: Path
    biogrid_tsv: Path
    dna_features: Path
    output_network: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [self.pombase_dir, self.biogrid_tsv, self.dna_features]:
            if not path.exists():
                raise ValueError(f"Required input path does not exist: {path}")
        self.output_network.parent.mkdir(parents=True, exist_ok=True)

    @property
    def ontology_dir(self) -> Path:
        return self.pombase_dir / "ontologies_and_associations"


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
def run(config: NetworkConfig) -> None:
    """Load GO gene2go, then collect network-level features for the coding-gene set."""
    coding_genes = read_coding_genes(config.dna_features)
    od = config.ontology_dir

    logger.info("Loading GO ontology data for network-level features")
    ontology_cfg = OntologyDataConfig(
        ontology_obo=od / "go-basic.obo",
        ontology_association_gaf=od / "gene_ontology_annotation.gaf.tsv",
        slim_terms_table=[od / "bp_go_slim_terms.tsv", od / "mf_go_slim_terms.tsv", od / "cc_go_slim_terms.tsv"],
    )
    _, _, _, gene2go, _, _ = load_ontology_data(
        ontology_cfg.load_data(),
        relationships={"is_a", "part_of"}, propagate_counts=True, load_obsolete=False, prt=None,
    )

    logger.info("Collecting network-level features")
    network_df = collect_network_level_features(config.pombase_dir, config.biogrid_tsv, gene2go, coding_genes)
    write_parquet(network_df, config.output_network)
    logger.success(f"Wrote {len(network_df)} gene rows to {config.output_network}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Collect network-level pombe features")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--biogrid-tsv", type=Path, required=True, help="BioGrid interaction table")
    parser.add_argument("--dna-features", type=Path, required=True, help="DNA-level features parquet (for coding-gene set)")
    parser.add_argument("--output", type=Path, required=True, dest="output_network", help="Output network-level features pickle")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run network-level collection, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = NetworkConfig(
            pombase_dir=args.pombase_dir, biogrid_tsv=args.biogrid_tsv,
            dna_features=args.dna_features, output_network=args.output_network,
        )
        config.validate()
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
