#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RNA-Level Feature Collection
==============================

Assembles mRNA abundance (Marguerat 2012) and mRNA kinetics (Harigaya 2016)
per coding gene. Reads the coding-gene set from the DNA-level pickle.

Input
-----
- A PomBase version directory (for gene metadata / id resolution)
- Literature tables (Marguerat 2012, Harigaya 2016)
- DNA-level features pickle (for the coding-gene set)

Output
------
- rna_features.pkl: per-gene RNA-level feature table (indexed by gene id)

Usage
-----
    python collect_rna_features.py \\
        --pombase-dir resources/external/pombase/2025-10-01 \\
        --literature-dir resources/literature \\
        --dna-features results/features/2025-10-01/_levels/dna_features.pkl \\
        --output results/features/2025-10-01/_levels/rna_features.pkl

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
from workflow.src.features.assembly import collect_rna_level_features, read_coding_genes


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class RnaConfig:
    """Inputs/outputs for RNA-level feature collection."""
    pombase_dir: Path
    literature_dir: Path
    dna_features: Path
    output_rna: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [self.pombase_dir, self.literature_dir, self.dna_features]:
            if not path.exists():
                raise ValueError(f"Required input path does not exist: {path}")
        self.output_rna.parent.mkdir(parents=True, exist_ok=True)

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
def run(config: RnaConfig) -> None:
    """Collect RNA-level features filtered to the DNA-level coding-gene set."""
    coding_genes = read_coding_genes(config.dna_features)
    logger.info("Collecting RNA-level features")
    rna_df = collect_rna_level_features(config.literature_dir, config.gene_meta_file, coding_genes)
    rna_df.to_pickle(config.output_rna)
    logger.success(f"Wrote {len(rna_df)} gene rows to {config.output_rna}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Collect RNA-level pombe features")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--literature-dir", type=Path, required=True, help="Directory of literature supplementary tables")
    parser.add_argument("--dna-features", type=Path, required=True, help="DNA-level features pickle (for coding-gene set)")
    parser.add_argument("--output", type=Path, required=True, dest="output_rna", help="Output RNA-level features pickle")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run RNA-level collection, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = RnaConfig(
            pombase_dir=args.pombase_dir, literature_dir=args.literature_dir,
            dna_features=args.dna_features, output_rna=args.output_rna,
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
