#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Phenotype-Level Feature Collection
====================================

Assembles FYPO viability, deletion-library essentiality, bar-seq fitness,
transposon insertion density, and CRISPRi growth phenotypes per coding gene.
Reads the coding-gene set from the DNA-level pickle.

NOTE: this level intentionally produces a DUPLICATE `DeletionLibrary_essentiality`
column (byte-faithful quirk carried through to the final matrix). The output is
pickled rather than TSV'd so the duplicate column name survives round-trip.

Input
-----
- A PomBase version directory (gene_viability.tsv)
- Curated deletion-library + essentiality-verification tables
- Literature tables (QianWenFeng/Koch, Guo 2013, Ishikawa 2024, Grech 2019)
- DNA-level features parquet (for the coding-gene set)

Output
------
- phenotype_features.parquet: per-gene phenotype feature table (indexed by gene id)

Usage
-----
    python collect_phenotype_features.py \\
        --pombase-dir resources/external/pombase/2025-10-01 \\
        --literature-dir resources/literature \\
        --deletion-library-xlsx resources/curated/deletion_library_categories.xlsx \\
        --essentiality-verification-csv resources/curated/essentiality_verification.csv \\
        --dna-features results/features/2025-10-01/_levels/dna_features.parquet \\
        --output results/features/2025-10-01/_levels/phenotype_features.parquet

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
    collect_phenotype_level_features,
    load_phyloP_and_divergence,
    read_coding_genes,
)


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PhenotypeConfig:
    """Inputs/outputs for phenotype-level feature collection."""
    pombase_dir: Path
    literature_dir: Path
    deletion_library_xlsx: Path
    essentiality_verification_csv: Path
    dna_features: Path
    output_phenotype: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [
            self.pombase_dir, self.literature_dir, self.deletion_library_xlsx,
            self.essentiality_verification_csv, self.dna_features,
        ]:
            if not path.exists():
                raise ValueError(f"Required input path does not exist: {path}")
        self.output_phenotype.parent.mkdir(parents=True, exist_ok=True)

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
def run(config: PhenotypeConfig) -> None:
    """Collect phenotype-level features filtered to the DNA-level coding-gene set."""
    coding_genes = read_coding_genes(config.dna_features)
    phyloP_and_divergence = load_phyloP_and_divergence(config.literature_dir, config.gene_meta_file)

    logger.info("Collecting phenotype-level features")
    phenotype_df = collect_phenotype_level_features(
        config.pombase_dir, config.deletion_library_xlsx, config.essentiality_verification_csv,
        config.literature_dir, config.gene_meta_file, coding_genes, phyloP_and_divergence,
    )
    write_parquet(phenotype_df, config.output_phenotype)
    logger.success(f"Wrote {len(phenotype_df)} gene rows to {config.output_phenotype}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Collect phenotype-level pombe features")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--literature-dir", type=Path, required=True, help="Directory of literature supplementary tables")
    parser.add_argument("--deletion-library-xlsx", type=Path, required=True, help="Curated deletion library categories xlsx")
    parser.add_argument("--essentiality-verification-csv", type=Path, required=True, help="Curated essentiality verification csv")
    parser.add_argument("--dna-features", type=Path, required=True, help="DNA-level features parquet (for coding-gene set)")
    parser.add_argument("--output", type=Path, required=True, dest="output_phenotype", help="Output phenotype-level features pickle")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run phenotype-level collection, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PhenotypeConfig(
            pombase_dir=args.pombase_dir, literature_dir=args.literature_dir,
            deletion_library_xlsx=args.deletion_library_xlsx,
            essentiality_verification_csv=args.essentiality_verification_csv,
            dna_features=args.dna_features, output_phenotype=args.output_phenotype,
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
