#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature Matrix Merge
======================

Outer-joins the six per-level feature pickles (DNA / RNA / protein /
evolutionary / network / phenotype) into the final per-coding-gene feature
matrix, matching the former monolithic collect_pombe_features.py output
byte-for-byte (including the intentional duplicate DeletionLibrary_essentiality
column, preserved via the parquet intermediates).

Input
-----
- The six per-level feature pickles
- A PomBase version directory (gene metadata for the gene_name column)

Output
------
- pombe_coding_gene_protein_features.tsv: one row per coding gene

Usage
-----
    python merge_features.py \\
        --pombase-dir resources/external/pombase/2025-10-01 \\
        --dna-features   results/features/2025-10-01/_levels/dna_features.parquet \\
        --rna-features   results/features/2025-10-01/_levels/rna_features.parquet \\
        --protein-features results/features/2025-10-01/_levels/protein_features.parquet \\
        --evolutionary-features results/features/2025-10-01/_levels/evolutionary_features.parquet \\
        --network-features results/features/2025-10-01/_levels/network_features.parquet \\
        --phenotype-features results/features/2025-10-01/_levels/phenotype_features.parquet \\
        --output results/features/2025-10-01/pombe_coding_gene_protein_features.tsv

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
from workflow.src.io import read_parquet, write_parquet
from workflow.src.features.assembly import load_gene_meta, merge_all_features


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class MergeConfig:
    """The six per-level pickles, gene metadata, and the final output path."""
    pombase_dir: Path
    dna_features: Path
    rna_features: Path
    protein_features: Path
    evolutionary_features: Path
    network_features: Path
    phenotype_features: Path
    output_features: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [
            self.pombase_dir, self.dna_features, self.rna_features, self.protein_features,
            self.evolutionary_features, self.network_features, self.phenotype_features,
        ]:
            if not path.exists():
                raise ValueError(f"Required input path does not exist: {path}")
        self.output_features.parent.mkdir(parents=True, exist_ok=True)

    @property
    def gene_meta_file(self) -> Path:
        """PomBase gene_IDs_names_products.tsv, used for the gene_name column."""
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
def run(config: MergeConfig) -> None:
    """Load the six per-level pickles, merge them, and write the final feature matrix."""
    gene_meta, _ = load_gene_meta(config.gene_meta_file)

    dna_df = read_parquet(config.dna_features)
    rna_df = read_parquet(config.rna_features)
    protein_df = read_parquet(config.protein_features)
    evolutionary_df = read_parquet(config.evolutionary_features)
    network_df = read_parquet(config.network_features)
    phenotype_df = read_parquet(config.phenotype_features)

    logger.info("Merging all feature groups")
    pombe_features = merge_all_features(dna_df, rna_df, protein_df, evolutionary_df, network_df, phenotype_df, gene_meta)

    pombe_features.to_csv(config.output_features, sep="\t", index=False)
    logger.success(f"Wrote {len(pombe_features)} gene records to {config.output_features}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Merge per-level pombe feature pickles into the final matrix")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--dna-features", type=Path, required=True, help="DNA-level features pickle")
    parser.add_argument("--rna-features", type=Path, required=True, help="RNA-level features pickle")
    parser.add_argument("--protein-features", type=Path, required=True, help="Protein-level features pickle")
    parser.add_argument("--evolutionary-features", type=Path, required=True, help="Evolutionary-level features pickle")
    parser.add_argument("--network-features", type=Path, required=True, help="Network-level features pickle")
    parser.add_argument("--phenotype-features", type=Path, required=True, help="Phenotype-level features pickle")
    parser.add_argument("--output", type=Path, required=True, dest="output_features", help="Output feature matrix path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the merge, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = MergeConfig(
            pombase_dir=args.pombase_dir,
            dna_features=args.dna_features,
            rna_features=args.rna_features,
            protein_features=args.protein_features,
            evolutionary_features=args.evolutionary_features,
            network_features=args.network_features,
            phenotype_features=args.phenotype_features,
            output_features=args.output_features,
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
