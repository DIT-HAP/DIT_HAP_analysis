#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein-Level Feature Collection
==================================

Assembles peptide-sequence features, protein abundance/turnover, AlphaFold
pLDDT statistics, and PFAM-domain counts per coding gene. Reads the
coding-gene set from the DNA-level pickle. This is the slowest level (it walks
the AlphaFold structure directory), so it is its own Snakemake rule.

Input
-----
- A PomBase version directory (peptide.fa, peptide_stats.tsv, gene metadata,
  RNA_metadata, Protein_features)
- An AlphaFold structure directory (.pdb.gz files)
- Literature tables (Christiano 2014)
- DNA-level features pickle (for the coding-gene set)

Output
------
- protein_features.pkl: per-gene protein-level feature table (indexed by gene id)

Usage
-----
    python collect_protein_features.py \\
        --pombase-dir resources/external/pombase/2025-10-01 \\
        --alphafold-dir /path/to/AlphaFold_Dataset \\
        --literature-dir resources/literature \\
        --dna-features results/features/2025-10-01/_levels/dna_features.pkl \\
        --output results/features/2025-10-01/_levels/protein_features.pkl

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
from workflow.src.features.assembly import collect_protein_level_features, load_gene_meta, read_coding_genes


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class ProteinConfig:
    """Inputs/outputs for protein-level feature collection."""
    pombase_dir: Path
    alphafold_dir: Path
    literature_dir: Path
    dna_features: Path
    output_protein: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [self.pombase_dir, self.alphafold_dir, self.literature_dir, self.dna_features]:
            if not path.exists():
                raise ValueError(f"Required input path does not exist: {path}")
        self.output_protein.parent.mkdir(parents=True, exist_ok=True)

    @property
    def gene_meta_file(self) -> Path:
        """PomBase gene_IDs_names_products.tsv, used for update_sysIDs()/uniprot2id."""
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
def run(config: ProteinConfig) -> None:
    """Collect protein-level features filtered to the DNA-level coding-gene set."""
    coding_genes = read_coding_genes(config.dna_features)
    _, uniprot2id = load_gene_meta(config.gene_meta_file)
    protein_meta = pd.read_csv(config.pombase_dir / "Protein_features" / "peptide_stats.tsv", sep="\t", index_col=0)

    logger.info("Collecting protein-level features")
    protein_df = collect_protein_level_features(
        config.pombase_dir, config.alphafold_dir, config.literature_dir,
        config.gene_meta_file, protein_meta, uniprot2id, coding_genes,
    )
    protein_df.to_pickle(config.output_protein)
    logger.success(f"Wrote {len(protein_df)} gene rows to {config.output_protein}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Collect protein-level pombe features")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--alphafold-dir", type=Path, required=True, help="AlphaFold structure directory (.pdb.gz files)")
    parser.add_argument("--literature-dir", type=Path, required=True, help="Directory of literature supplementary tables")
    parser.add_argument("--dna-features", type=Path, required=True, help="DNA-level features pickle (for coding-gene set)")
    parser.add_argument("--output", type=Path, required=True, dest="output_protein", help="Output protein-level features pickle")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run protein-level collection, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = ProteinConfig(
            pombase_dir=args.pombase_dir, alphafold_dir=args.alphafold_dir, literature_dir=args.literature_dir,
            dna_features=args.dna_features, output_protein=args.output_protein,
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
