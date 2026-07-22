#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Evolutionary-Level Feature Collection
=======================================

Assembles ortholog/paralog counts, evolutionary rate, and phyloP/divergence
scores per coding gene. Reads the coding-gene set from the DNA-level pickle.

Input
-----
- A PomBase version directory (curated_orthologs, gene metadata)
- An Ensembl paralog export TSV
- Literature tables (Rhind 2011, Grech 2019)
- DNA-level features parquet (for the coding-gene set)

Output
------
- evolutionary_features.parquet: per-gene evolutionary feature table (indexed by gene id)

Usage
-----
    python collect_evolutionary_features.py \\
        --pombase-dir resources/external/pombase/2025-10-01 \\
        --literature-dir resources/literature \\
        --ensembl-paralogs-tsv resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv \\
        --dna-features results/features/2025-10-01/_levels/dna_features.parquet \\
        --output results/features/2025-10-01/_levels/evolutionary_features.parquet

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
from workflow.src.features.assembly import (
    collect_evolutionary_level_features,
    load_phyloP_and_divergence,
    read_coding_genes,
)


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class EvolutionaryConfig:
    """Inputs/outputs for evolutionary-level feature collection."""
    pombase_dir: Path
    literature_dir: Path
    ensembl_paralogs_tsv: Path
    dna_features: Path
    output_evolutionary: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [self.pombase_dir, self.literature_dir, self.ensembl_paralogs_tsv, self.dna_features]:
            if not path.exists():
                raise ValueError(f"Required input path does not exist: {path}")
        self.output_evolutionary.parent.mkdir(parents=True, exist_ok=True)

    @property
    def gene_meta_file(self) -> Path:
        """PomBase gene_IDs_names_products.tsv, used for update_sysIDs()."""
        return self.pombase_dir / "Gene_metadata" / "gene_IDs_names_products.tsv"


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
def run(config: EvolutionaryConfig) -> None:
    """Collect evolutionary-level features filtered to the DNA-level coding-gene set."""
    coding_genes = read_coding_genes(config.dna_features)
    phyloP_and_divergence = load_phyloP_and_divergence(config.literature_dir, config.gene_meta_file)

    logger.info("Collecting evolutionary-level features")
    evolutionary_df = collect_evolutionary_level_features(
        config.pombase_dir, config.ensembl_paralogs_tsv, config.literature_dir,
        config.gene_meta_file, coding_genes, phyloP_and_divergence,
    )
    write_parquet(evolutionary_df, config.output_evolutionary)
    logger.success(f"Wrote {len(evolutionary_df)} gene rows to {config.output_evolutionary}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Collect evolutionary-level pombe features")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--literature-dir", type=Path, required=True, help="Directory of literature supplementary tables")
    parser.add_argument("--ensembl-paralogs-tsv", type=Path, required=True, help="Ensembl paralog export table")
    parser.add_argument("--dna-features", type=Path, required=True, help="DNA-level features parquet (for coding-gene set)")
    parser.add_argument("--output", type=Path, required=True, dest="output_evolutionary", help="Output evolutionary-level features pickle")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run evolutionary-level collection, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = EvolutionaryConfig(
            pombase_dir=args.pombase_dir, literature_dir=args.literature_dir,
            ensembl_paralogs_tsv=args.ensembl_paralogs_tsv,
            dna_features=args.dna_features, output_evolutionary=args.output_evolutionary,
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
